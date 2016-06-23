#include "fmm.h"

#ifdef GRAVITY_FMM

static void prepare_fmm();

static int find_leaf_size(const int, const int, const int);
static int find_next_level(const int, int);
static void copy_leafs(const int, const int *, const int *)

uint64_t NNodes = 0, Max_Nodes = 0;
extern struct FMM_Node FMM = { NULL };
static struct FMM_Node fmm = { NULL };
#pragma omp threadprivate(fmm)
static omp_lock_t Node_Lock;
static Leaf_Nodes[] = NULL;

/*
 * Build the FMM tree in parallel. 
 * This is an exercise in OpenMP tasking as it is MPI task local.
 * At this point we have NTop_Nodes containing npart starting at First_Part.
 * We find the leafs of the top node tree from the peano key and then 
 * kick-off the tree build and the FMM P2M on the go.
 *
 */

void Gravity_FMM_Build_P2M()
{
	Profile("Grav FMM Build & P2M");

	#pragma omp single
	NNodes = 0;
	
// for (;;)
		#pragma omp single
		prepare_fmm();

		//#pragma omp taskloop
		//for (int i = 0; i < NTop_Nodes; i++) {
		
			#pragma omp taskyield

	//		if (D[i].TNode.Target < 0)  // not local
	//			continue;
	
			struct FMM_Node fmm = reserve_memory();
			int nNodes = 0;

			/* Works because buffer has space of Npart_Total_Max size_t's */

			size_t nBytes = Task.Npart_Total_Max*sizeof(int);

			int *leaf2part = Get_Thread_Safe_Buffer(nBytes); 
			int *leaf2node = leaf2part + 0.5 * nBytes; 
	
			const int first_part = D[i].TNode.First_Part;
			const int last_part = first_part + D[i].TNode.Npart - 1;
			const int tnode_lvl = D[i].TNode.Level + 1;

			int ipart = find_leaf_size(first_part, last_part, tnode_lvl);;
			int lvl = find_next_level(ipart, tnode_lvl);
	
			leaf2part[0] = first_part;
			leaf2node[0] = build_tree(first_part, ipart, lvl, fmm, nNodes);

			int nLeafs = 1;

			while (ipart <= last_part) { // all particles in the top node
	
				int npart = find_leaf_size(ipart, last_part, lvl);

				leaf2part[nLeafs] = ipart;
				leaf2node[nLeafs] = build_tree(ipart, npart, lvl, fmm, nNodes);

				nLeafs++;

				if (nLeafs > 1) {

					#pragma omp task untied
					fmm_p2m(leaf2node[nLeafs-1], jpart, fmm);
				}

				ipart += npart;

				lvl = find_next_level(ipart, D[i].TNode.Level + 1);

			} // while

			#pragma omp task untied
			fmm_p2m(leaf2node[nLeafs], jpart, fmm);
	
			#pragma omp critical
			D[i].TNode.First_Leaf = copy_leafs(i,nLeafs, leaf2part, leaf2node);

		} // for i

		#pragma omp barrier

	} // forever

	Profile("Grav FMM Build");

	return ;
}

/*
 * Find leaf vector of particle ipart at level lvl and return the number of
 * particles in the leave.
 */

static int find_leaf_size(const int ipart, const int last, const int lvl)
{
	if (ipart == last) // last particle -> single leaf
		return last + 1;

	int jmax = MIN(last, ipart + VECTOR_SIZE) + 1;

	int npart = 0; 

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

		if (npart <= VECTOR_SIZE) // leaf found
			break;

		mask <<= 3;
			
		lvl++;

		//Assert(lvl < 30, "%d %d %d", lvl, ipart, i);
		
	} // forever

	return npart; // next leaf
}

static int find_next_level(const int ipart, int lvl)
{
	mask = ((peanoKey) 0x7) << (3*lvl);

	while ((P.Key[ipart] & mask) == (P.Key[ipart-1] & mask)) {

		lvl++;
		mask <<= 3;
	}

	return lvl;
}

static void copy_leafs(const int n, const int *leafs, const int *nodes)
{
	int dest = 0;

	dest = NLeafs;
	NLeafs += n;

	memcpy(&Leaf2Part[dest], leafs, n*sizeof(*leafs));
	memcpy(&Leaf2Node[dest], nodes, n*sizeof(*nodes));

	return dest;
}

static void fmm_tree_build()
{

fmm = FMM;

	int first_vec = D[0].TNode.First_Vec;
	int top_level = D[0].TNode.Level;

	peanoKey last_key = create_first_node(Vec[first_vec], 0, D[0].TNode.Level,
									      first_vec);
	NNodes = 1;

	int last_parent = 0;

	for (int ivec = first_vec+1; ivec < D[0].TNode.NVec; ivec++) {
	
		int ipart = Vec[ivec];

		peanoKey key = P.Key[ipart];

		key >>= 3 * top_level;

		int node = 0;        // running node
		int lvl = top_level; // counts current level
		int parent = node;   // parent of current node

		while (lvl < N_PEANO_TRIPLETS) {

			if (particle_is_inside_node(key, lvl, node)) { // open node	

				if (fmm[node].Npart == 1) { // refine
				
					int jpart = -fmm[node].DNext -1;

					fmm[node].DNext = 0;

					int new_node = nNodes; // is a son of "node"

					create_node_from_particle(jpart, node, last_key, lvl+1,
											  new_node, ivec);
					nNodes++;

					last_key >>= 3;
				}

				parent = node;

				node++; // decline
				lvl++;
				key >>= 3;

			} else { // skip right to next node
	
				if (fmm[node].DNext == 0 || node == nNodes - 1)
					break; // reached end of branch

				node += fmax(1, tree[node].DNext);
			}
		} // while lvl

		if (fmm[node].DNext == 0)	// set DNext for internal node
			fmm[node].DNext = nNodes - node;	// only delta

		int new_node = nNodes; // is a sibling of "node"

		create_node_from_particle(ipart, parent, key, lvl, new_node, ivec);

		nNodes++;

		last_key = key >> 3;

		last_parent = parent;
			
	} // for i


}

static void prepare_fmm()
{
	if (NNodes != 0) { // tree build aborted, increase mem

		Max_Nodes = TREE_ENLARGEMENT_FACTOR * Max_Nodes;
		
		printf("(%d:%d) Increased tree memory to %6.1f MB, "
			"max %10d nodes, ratio %4g \n", Task.Rank, Task.Thread_ID, 
			Max_Nodes * sizeof(*Tree)/1024.0/1024.0, Max_Nodes, 
			(double) Max_Nodes/Task.Npart_Total); 
	}

	Alloc_Gravity_FMM(Max_Nodes, FMM);
		
	NNodes = 0;

	return ;
}

/*
 * The first node in the subtree needs special treatment. Its properties are 
 * identical with the topnode from the domain decomposition.
 */

static peanoKey create_first_node(const int first_part,
		const int tnode_idx, const int top_level, const int ivec)
{
	peanoKey key = P.Key[first_part];

	key >>= 3 * top_level;

	create_node_from_particle(first_part, 0, key, top_level, 0, ivec);

	fmm[0].Pos[0] = D[tnode_idx].TNode.Pos[0]; // get top node position 
	fmm[0].Pos[1] = D[tnode_idx].TNode.Pos[1]; // because parent node 
	fmm[0].Pos[2] = D[tnode_idx].TNode.Pos[2]; // did not exist

	fmm[0].DUp = tnode_idx; // correct up pointer to lead to topnode

	node_set(TOP, 0);

	return key >> 3;
}


/*
 * For particle and node to overlap the peano key triplet at this tree level 
 * has to be equal. Hence the tree cannot be deeper than the PH key
 * resolution, which for our 128 bits length is 42. This corresponds to 
 * distances less than 2^-42, small enough for single precision.
 */

static inline bool particle_is_inside_node(const peanoKey key, const int lvl,
										   const int node)
{
	int part_triplet = key & 0x7;

	int node_triplet = key_fragment(node);

	return node_triplet == part_triplet;
}

/*
 * We always add nodes at the end of the subtree. Particles are named 
 * negative and offset by one to leave DNext=0 indicating unset. We assume
 * the peano key has reversed triplet order and the least significant 3 bits 
 * carry the triplet at level "lvl".
 */

static inline void create_node_from_particle(const int ipart, const int parent,
											 const peanoKey key, const int lvl,
											 const int node, const int ivec)
{
	fmm[node].DNext = -ivec - 1;

	int keyfragment = (key & 0x7) << 6;

	fmm[node].Bitfield = lvl | keyfragment | (1UL << 9);

	const int sign[3] = { -1 + 2 * (P.Pos[0][ipart] > tree[parent].Pos[0]),
					      -1 + 2 * (P.Pos[1][ipart] > tree[parent].Pos[1]),
		                  -1 + 2 * (P.Pos[2][ipart] > tree[parent].Pos[2]) };

	Float size = Domain.Size / (1 << lvl);

	fmm[node].Pos[0] = fmm[parent].Pos[0] + sign[0] * size * 0.5;
	fmm[node].Pos[1] = fmm[parent].Pos[1] + sign[1] * size * 0.5;
	fmm[node].Pos[2] = fmm[parent].Pos[2] + sign[2] * size * 0.5;

	fmm[node].DUp = node - parent;

	P.Tree_Parent[ipart] = node;

	return ;
}

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

static inline void node_set(const enum Tree_Bitfield bit, const int node)
{
	tree[node].Bitfield |= 1UL << bit;

	return ;
}



#endif // GRAVITY_FMM
