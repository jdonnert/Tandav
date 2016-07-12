#include "fmm.h"

#ifdef GRAVITY_FMM

#define NODES_PER_PARTICLE 0.6
#define TREE_ENLARGEMENT_FACTOR 1.2

#define DEBUG_FMM

static void prepare_fmm();
static struct FMM_Node reserve_fmm_memory(const int, int *);

static int find_leaf_size(const int, const int, const int);
static int find_next_level(const int, int);
static int copy_leafs(const int, const int *, const int *);
static void set_node_position(const int, const int, const int, const int,
							  const int, struct FMM_Node fmm);
static int fmm_build_branch(const int, const int, const int, const int,
							struct FMM_Node, int);
static void fmm_p2m(const int, int, const int, struct FMM_Node);
static bool particle_is_inside_node(const peanoKey, const int);

void test_tree (struct FMM_Node f, int *leaf2part, int *leaf2node,
				const int nLeafs, const int nNodes, const int tnode);

static struct FMM_Node fmm = { NULL }; // we build the tree in here
#pragma omp threadprivate(fmm)

/*
 * Build the FMM tree. This is an exercise in OpenMP tasking as it is MPI task
 * local. We safe to leaf2X arrays in the thread safe buffer and copy it back.
 * At this point we have NTop_Nodes containing npart starting at First_Part.
 * We find the leafs of the top node tree from the peano key and then
 * kick-off the build of a branch of the tree and the FMM P2M.
 *
 */

void Gravity_FMM_Build()
{
	Profile("Grav FMM Build");

	#pragma omp single
	NNodes = 0;

	size_t nBytes = 2 * Task.Npart_Total_Max*sizeof(int);
	int *leaf2part = Get_Thread_Safe_Buffer(nBytes);
	int *leaf2node = leaf2part + (nBytes >> 1);

 	do {

		#pragma omp barrier
		
		#pragma omp single
		prepare_fmm();

		#pragma omp single // for
		for (int i = 0; i < NTop_Nodes; i++) {
			
printf("Tnode level = %d first_part =%d %d\n", D[i].TNode.Level, D[i].TNode.First_Part, NTop_Nodes);

		//	if (D[i].TNode.Target != -Task.Rank-1) // on other task
		//		continue;

			int nNodes = 0;
			int offset = 0;

			struct FMM_Node fmm = reserve_fmm_memory(i, &offset); // FMM[offset]

			if (NNodes >= Max_Nodes) // no more memory, force rebuild
				continue;

			int nLeafs = 0;

			const int last_part = D[i].TNode.First_Part + D[i].TNode.Npart-1;
			const int min_lvl = D[i].TNode.Level + 1;

			int ipart = D[i].TNode.First_Part;
			int level = min_lvl;

			while (ipart <= last_part) { // all particles in the top node

				int npart = find_leaf_size(ipart, last_part, level);

				level = find_next_level(ipart+npart, min_lvl); // npart>0

				nNodes = fmm_build_branch(ipart, i, npart, level,
										  fmm, nNodes); // leaf now at end

				int node = nNodes - 1;

				fmm.Leaf_Ptr[node] = -nLeafs - 1; // insert leaf

				leaf2node[nLeafs] = node + offset; // correct offset
				leaf2part[nLeafs] = ipart;

//				#pragma omp task untied
				fmm_p2m(ipart, nNodes+offset, npart, fmm);

				nLeafs++;

				ipart += npart;

			} // while ipart

			fmm.DUp[0] = i; // points to topnode

			/* from here on D.TNode.First_Leaf union points to leafs */

printf("Tnode level = %d first_part =%d \n", D[i].TNode.Level, D[i].TNode.First_Part);
			D[i].TNode.Target = offset;
			D[i].TNode.First_Leaf = copy_leafs(nLeafs, leaf2part, leaf2node);
			D[i].TNode.NLeafs = nLeafs;
			D[i].TNode.Mass = fmm.Mass[0];
			D[i].TNode.CoM[0] = fmm.CoM[0][0];
			D[i].TNode.CoM[1] = fmm.CoM[1][0];
			D[i].TNode.CoM[2] = fmm.CoM[2][0];
			D[i].TNode.Dp[0] = fmm.Dp[0][0];
			D[i].TNode.Dp[1] = fmm.Dp[1][0];
			D[i].TNode.Dp[2] = fmm.Dp[2][0];
			
			test_tree(fmm, leaf2part, leaf2node, nLeafs, nNodes,i); //DEBUG_FMM
		} // for i

	} while (NNodes >= Max_Nodes); // while tree mem too small
exit(0);

	Profile("Grav FMM Build");

	return ;
}

static void prepare_fmm()
{
	if (NNodes > 0) { // tree build aborted, increase memory

		Max_Nodes *= TREE_ENLARGEMENT_FACTOR;

		printf("(%d:%d) Increased FMM tree memory to %6.1f MB, "
			"max %10d nodes, ratio %4g \n", Task.Rank, Task.Thread_ID,
			Max_Nodes * sizeof_FMM/1024.0/1024.0, Max_Nodes,
			(double) Max_Nodes/Task.Npart_Total);

		for (int i = 0; i < NTop_Nodes; i++) { // restore top nodes
			
			if (D[i].TNode.Target < 0) // untouched or other rank
				continue;
			
			D[i].TNode.Target = -Task.Rank - 1;

			int first_leaf = D[i].TNode.First_Leaf;
			D[i].TNode.First_Part = Leaf2Part[first_leaf];

		}
			memset(Leaf2Part, 0, Task.Npart_Total * sizeof(*Leaf2Part));
			memset(Leaf2Node, 0, Task.Npart_Total * sizeof(*Leaf2Node));
	}

	Gravity_FMM_Free(FMM);

	FMM = Alloc_FMM_Nodes(Max_Nodes);


	Print_Memory_Usage();

	NNodes = 0;

	return ;
}

/*
 * This returns a collection of pointers to the global node structure FMM,
 * starting at the first node for this subtree. The offset of the pointers
 * to the global FMM structure are written into *offset
 */

static struct FMM_Node reserve_fmm_memory(const int tnode, int *offset)
{
	int nReserved = ceil(D[tnode].TNode.Npart * NODES_PER_PARTICLE);

	omp_set_lock(&NNodes_Lock);

	int i = NNodes;

	NNodes += nReserved;

	omp_unset_lock(&NNodes_Lock);

	struct FMM_Node f = { NULL };

	if (NNodes > Max_Nodes) // abort tree build and alloc more memory
		return f;
	
	f.DNext = &FMM.DNext[i];
	f.Bitfield = &FMM.Bitfield[i];
	f.DUp = &FMM.DUp[i];
	f.Npart = &FMM.Npart[i];
	f.Mass = &FMM.Mass[i];

	for (int j = 0; j < 3; j++) {

		f.CoM[j] = &FMM.CoM[j][i];
		f.Dp[j] = &FMM.Dp[j][i];

#ifdef FMM_SAVE_NODE_POS
		f.Pos[j] = &FMM.Pos[j][i];
#endif
	}

	*offset = i;

	return f;
}

/*
 * Find leaf vector of particle ipart at level lvl and return the number of
 * particles in the leaf. At a leaf, not more than VECTOR_SIZE particles have
 * the same key at all levels above the level of the leaf.
 * If bunch of particles are closer than the PH key precision, we increase
 * the size of the leaf until the level is fits into the leaf.
 */

static int find_leaf_size(const int ipart, const int last_part,
						  const int min_lvl)
{
	if (ipart > last_part - VECTOR_SIZE) // manage end
		return last_part - ipart + 1;

	const int max_lvl = N_PEANO_TRIPLETS-1;

	int npart = 0;
	int lvl = 0;
	int size = VECTOR_SIZE;

	do {

		npart = 0;
		lvl = min_lvl;

		int jmax = MIN(last_part, ipart + size) + 1;

		npart = 0;

		peanoKey mask = ((peanoKey) 0x7) << (3*lvl);

		for (;;) { // loop over level

			npart = 1; // account for ipart

			const peanoKey triplet_i = P.Key[ipart] & mask;

			for (int jpart = ipart + 1; jpart < jmax; jpart++) {

				peanoKey triplet_j = P.Key[jpart] & mask;

				if (triplet_j != triplet_i)
					break;

				npart++;
			}

			if (npart <= size || lvl == max_lvl)
				break; // leaf found: npart with same key at lvl

			mask <<= 3;

			lvl++;

		} // forever

		size++; // create leaf larger than VECTOR_SIZE

	} while (lvl == max_lvl); // repeat until we fit

	Assert(npart < 1024, "I found a very large leaf (npart=%d) at ipart=%d. "
				"That means your simulation has a very large spatial "
				"dynamic range. Many particle distances are only a factor "
				"of 8 away from the single precision limit and the "
				"algorithm becomes inaccurate ! Use DOUBLE_PRECISION to "
				"speed up the computation and fix the issue", npart, ipart);

	return npart; // until first particle of next leaf
}

static int find_next_level(const int ipart, int lvl)
{
	peanoKey mask = ((peanoKey) 0x7) << (3*lvl);

	while ((P.Key[ipart] & mask) == (P.Key[ipart-1] & mask)) {

		lvl++;
		mask <<= 3;
	}

	return lvl;
}

/*
 * Move the leaf to particle and node pointers into global memory.
 */

static int copy_leafs(const int n, const int *leafs, const int *nodes)
{
	int dest = 0;

	omp_set_lock(&NLeafs_Lock);

	dest = NLeafs;
	NLeafs += n;

	omp_unset_lock(&NLeafs_Lock);

	memcpy(&Leaf2Part[dest], leafs, n*sizeof(*leafs));
	memcpy(&Leaf2Node[dest], nodes, n*sizeof(*nodes));

	return dest;
}

/*
 * Build a branch of the FMM tree in *fmm.
 * As we know the level of the tree leaf, we refine and create nodes until we
 * reach that level.
 * We always add nodes at the end of the subtree. In a leaf, the
 * DNext/Leaf_Ptr union holds the index to leaf2particle. In an internal
 * node it holds the number of nodes to the next node.
 * The leaf index is stored negative and offset by one to leave 0 indicating
 * unset. We assume the peano key has reversed triplet order and the least
 * significant 3 bits carry the triplet at level "lvl".
 */

static int fmm_build_branch(const int ipart, const int tnode,
		const int npart, const int leaf_lvl, struct FMM_Node fmm, int nNodes)
{
	int top_level = D[tnode].TNode.Level;
	peanoKey key = P.Key[ipart] >> (3*top_level);

	int lvl = top_level; // counts current level

	int node = 0;        // running node
	int parent = 0;

	while (lvl < leaf_lvl+1) { // build the complete branch

		if (node == nNodes) { // end of branch, create node & decline

			fmm.DUp[node] = node - parent;
			fmm.Bitfield[node] = (lvl << 3) | ((key & 0x7)) ; // lvl < 2^6
			fmm.Bitfield[node] |= (int)(node == 0) << 9; // mark topnode

			nNodes++;

			set_node_position(ipart, node, tnode, parent, lvl, fmm);

			lvl++;
			key >>= 3;
			parent = node;

		} else if (particle_is_inside_node(key, fmm.Bitfield[node])) {

			lvl++; // decline
			key >>= 3;
			parent = node;

		} else if (fmm.DNext[node] == 0) { // open new branch

			fmm.DNext[node] = nNodes - node;
		}

		node += imax(1, fmm.DNext[node]); // walk downward or forward
	}


	return nNodes;
}

/*
 * For particle and node to overlap the peano key triplet down to this tree
 * level has to be equal. Hence the tree cannot be deeper than the PH key
 * resolution, which for our 128 bits length is 42. This corresponds to
 * distances less than 2^-42, small enough for single precision. Double
 * precision positions would require higher precision.
 */

static bool particle_is_inside_node(const peanoKey key, const int bitfield)
{
	int part_triplet = key & 0x7;

	int node_triplet = bitfield & 0x7;

	return node_triplet == part_triplet;
}



static void set_node_position(const int ipart, const int node,
		const int tnode, const int parent, const int lvl, struct FMM_Node fmm)
{
#ifdef FMM_SAVE_NODE_POS

	if (node == 0) {

		fmm.Pos[0][0] = D[tnode].TNode.Pos[0]; // get top node position
		fmm.Pos[1][0] = D[tnode].TNode.Pos[1]; // because parent node
		fmm.Pos[2][0] = D[tnode].TNode.Pos[2]; // does not exist

	} else {

		int sign[3] = { -1 + 2 * (int)(P.Pos[0][ipart] > fmm.Pos[0][parent]),
				        -1 + 2 * (int)(P.Pos[1][ipart] > fmm.Pos[1][parent]),
		    	        -1 + 2 * (int)(P.Pos[2][ipart] > fmm.Pos[2][parent])};

		Float half_size = Domain.Size / (1 << (lvl+1));

		fmm.Pos[0][node] = fmm.Pos[0][parent] + sign[0] * half_size;
		fmm.Pos[1][node] = fmm.Pos[1][parent] + sign[1] * half_size;
		fmm.Pos[2][node] = fmm.Pos[2][parent] + sign[2] * half_size;
	}
#endif

	return ;
}


/*
 * Particle to Multipole sweep.
 * At the leaf node, the multipole and various other quantities in
 * the leaf are constructed. Then the fmm branch is walked upwards from the
 * leaf, adding the leaf values to the parent nodes. This is an O(n) operation
 * that vectorizes and can use fused multiply-add operations.
 * (Dehnen 2002, sect 3.1, Yokota 2012).
 */

static void fmm_p2m(const int beg, int node, const int npart,
					struct FMM_Node fmm)
{
	const int end = beg + npart;

	Float leaf_mass = 0,
		  leaf_CoM[3] = { 0 },
		  leaf_Dp[3] = { 0 };

	for (int ipart = beg; ipart < end; ipart++) { // vectorizes

		leaf_mass += P.Mass[ipart];

		leaf_CoM[0] += P.Mass[ipart] * P.Pos[0][ipart]; // fmadd
		leaf_CoM[1] += P.Mass[ipart] * P.Pos[1][ipart];
		leaf_CoM[2] += P.Mass[ipart] * P.Pos[2][ipart];

		leaf_Dp[0] += P.Mass[ipart] * P.Vel[0][ipart];
		leaf_Dp[1] += P.Mass[ipart] * P.Vel[1][ipart];
		leaf_Dp[2] += P.Mass[ipart] * P.Vel[2][ipart];
	}

	leaf_CoM[0] /= leaf_mass;
	leaf_CoM[1] /= leaf_mass;
	leaf_CoM[2] /= leaf_mass;

	int lvl = fmm.Bitfield[node] >> 3;

	for (int i = lvl; i >= 0; i--) {

		#pragma omp atomic update
		fmm.Npart[node] += npart;

		#pragma omp atomic update
		fmm.Mass[node] += leaf_mass;

		#pragma omp atomic update
		fmm.CoM[0][node] += leaf_CoM[0];
		#pragma omp atomic update
		fmm.CoM[1][node] += leaf_CoM[1];
		#pragma omp atomic update
		fmm.CoM[2][node] += leaf_CoM[2];

		#pragma omp atomic update
		fmm.Dp[0][node] += leaf_Dp[0];
		#pragma omp atomic update
		fmm.Dp[1][node] += leaf_Dp[1];
		#pragma omp atomic update
		fmm.Dp[2][node] += leaf_Dp[2];

		node -= fmm.DUp[node];
	} // for lvl

	return ;
}


void test_tree (struct FMM_Node f, int *leaf2part, int *leaf2node,
				const int nLeafs, const int nNodes, const int tnode)
{
#ifdef DEBUG_FMM
	printf("testing tree for consistency ... \n");

	int first_part = leaf2part[0];
	
	/*for (int i = first_part; i < first_part+20; i++)
		Print_Int_Bits(P.Key[i], 64, 0);

	for (int k =0;  k < nNodes; k++) {

		printf("%0d DNext=%d DUp=%d npart=%d lvl=%d key=%d%d%d ",
			k, f.DNext[k], f.DUp[k],f.Npart[k],
			(f.Bitfield[k] & 63 << 3) >> 3, (f.Bitfield[k] & (0x4)) >>2,
			(f.Bitfield[k] & (0x2)) >>1, f.Bitfield[k] & 0x1 );

		if (f.DNext[k] < 0) {

			int ll = -f.Leaf_Ptr[k] - 1;

			printf("leaf=%d node=%d start_part=%d \n", ll,
					leaf2node[ll], leaf2part[ll]);
		} else
			printf("\n");
	}*/

	int is_topnode = ((f.Bitfield[0] >> 9) & (1));

	Assert(is_topnode == 1, "Node 0 isn't marked as topnode");
	Assert(tnode == f.DUp[0], "TNode idx seems wrong");

	Assert(D[tnode].TNode.NLeafs == nLeafs, "nLeafs in tnode wrong %d - %d"
			, D[tnode].TNode.NLeafs, nLeafs);

	Assert(D[tnode].TNode.Mass == f.Mass[0], "Topnode mass wrong node=%d",0);
	Assert(D[tnode].TNode.CoM[0] == f.CoM[0][0], "CoM0 wrong %d", 0);
	Assert(D[tnode].TNode.CoM[1] == f.CoM[1][0], "CoM1 wrong %d", 0);
	Assert(D[tnode].TNode.CoM[2] == f.CoM[2][0], "CoM2 wrong %d", 0);

	for (int node = 1; node < nNodes; node++) {
		
		int lvl = (f.Bitfield[node] >> 3) & 63;
		int keyfrag = f.Bitfield[node] & 0x7;

		if (f.DNext[node] >= 0) { // internal node

			int dad = node - f.DUp[node];
			int dad_lvl = (f.Bitfield[dad] >> 3) & 63;

			Assert(dad_lvl+1 == lvl, "dad lvl wrong, node=%d, dad=%d"
					" dad_lvl = %d lvl=%d", node, dad,dad_lvl, lvl);

			int sister = node + f.DNext[node];

			if (sister != node) {

				int sister_lvl = (f.Bitfield[sister] >> 3) & 63;
				Assert(lvl == sister_lvl, "sister_lvl wrong, node=%d", node);


				int sister_keyfrag = f.Bitfield[sister] & 0x7;
				Assert(keyfrag < sister_keyfrag,"Sisters Keyfrag wrong %d",
						node);
			}

		} else if (f.DNext < 0) { // particle bundle

			int leaf = -f.DNext[node] - 1;

			Assert(node == leaf2node[leaf], "leaf2node or leaf ptr wrong, %d",
					node);

			int ipart = leaf2part[leaf];

			peanoKey keyfragment = f.Bitfield[node] & 0x7;
			peanoKey p_key = (P.Key[ipart] >> (3*lvl)) & 0x7;

			Assert(keyfragment == p_key, "key triplets wrong, %d", node);

			int jpart = leaf2part[leaf+1];
			int npart = jpart-ipart;

			Assert(f.Npart[node] == npart, "npart wrong, %d", node);
		}
	}

	printf("done \n");
#endif
	return ;
}

#endif // GRAVITY_FMM

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
