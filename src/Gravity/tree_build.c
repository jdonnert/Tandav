#include "../globals.h"
#include "gravity.h"
#include "../domain.h"
#include "../peano.h"
#include "../timestep.h"

#ifdef GRAVITY_TREE

#define NODES_PER_PARTICLE 0.6

static void gravity_tree_init();
static void transform_bunch_into_top_node(const int, int*, int*);
static int reserve_tree_memory(const int, const int);
static int build_subtree(const int, const int, const int);
static int finalise_subtree(const int, const int, int );
static inline bool particle_is_inside_node(const peanoKey,const int,const int);
static inline void add_particle_to_node(const int, const int);
static inline void create_node_from_particle(const int, const int,
											 const peanoKey, const int,
											 const int);
static peanoKey create_first_node(const int, const int, const int);
static void collapse_last_branch(const int, const int, const int, int*);
static inline int key_fragment(const int);
static inline void node_set(const enum Tree_Bitfield, const int);

int NNodes = 0;
static int Max_Nodes = 1024;

struct Tree_Node  *Tree = NULL; // pointer to all nodes

static struct Tree_Node *tree = NULL; //  build pointer "*Tree" or "*Buffer"
#pragma omp threadprivate(tree)

/*
 * This builds the tree in parallel, particles are assumed PH ordered. 
 * First every local bunch is converted into a topnode. Then we build the 
 * corresponding tree either in the openmp buffer, or directly inside the
 * Tree memory. 
 * A subtree is build starting from the top node in target pointer "*tree". 
 * Every subtree allocated a fixed amount of memory. 
 * If the number of particles in that top node is <= 8, the subtree is 
 * discarded and the topnode target points directly to the particles and 
 * the tree walk will use particles directly from the topnode. 
 * At last the top node is broadcasted.
 */

void Gravity_Tree_Build()
{
	Profile("Build Gravity Tree");

	#pragma omp single
	{

	NTop_Nodes = NBunches;

	NNodes = 0;

	Tree = Realloc(Tree, Max_Nodes * sizeof(*Tree), "Tree");

	} // omp single

	const size_t buf_thres = Task.Buffer_Size/sizeof(*Tree);

	#pragma omp for schedule(static,1)
	for (int i = 0; i < NTop_Nodes; i++) {

		int src = D[i].Bunch.Target;

		if (D[i].Bunch.Is_Local) {

			src = Task.Rank;

			int first_part = 0, level = 0;

			transform_bunch_into_top_node(i, &level, &first_part);

			bool build_in_buffer = D[i].TNode.Npart < buf_thres;

			if (build_in_buffer) {

				tree = Get_Thread_Safe_Buffer(Task.Buffer_Size);

			} else { // build in *Tree

				int nReserved = ceil(D[i].TNode.Npart * NODES_PER_PARTICLE);

				#pragma omp critical
				D[i].TNode.Target = reserve_tree_memory(i, nReserved);

				tree = &Tree[D[i].TNode.Target];
			}

			int nNeeded = build_subtree(first_part, i, level);

			if (build_in_buffer) { // copy buffer to Tree, clear buffer

				#pragma omp critical
				D[i].TNode.Target = reserve_tree_memory(i, nNeeded);

				size_t nBytes = nNeeded * sizeof(*Tree);

				memcpy(&Tree[D[i].TNode.Target], tree, nBytes);
			}

			int last_part = first_part + D[i].TNode.Npart;

			if (D[i].TNode.Target > 0) // correct particle parent pointer
				for (int ipart = first_part; ipart < last_part; ipart++)
					P[ipart].Tree_Parent += D[i].TNode.Target;
			else if (D[i].TNode.Target < 0)
				for (int ipart = first_part; ipart < last_part; ipart++)
					P[ipart].Tree_Parent = -i - 1; // top node w/o tree

		} // if Bunch local

	/*	MPI_Request *request = NULL;
		
		float *target = &D[i].TNode.Pos[0];
		int nBytes = (&D[i].TNode.Dp[2] - target) + sizeof(float);

		MPI_Ibcast(target, nBytes, MPI_BYTE, src, MPI_COMM_WORLD, request); */
	}

	rprintf("Tree build: %d of %d Nodes used (%g MB)\n",
			NNodes, Max_Nodes, Max_Nodes*sizeof(*Tree)/1024.0/1024);

#ifdef DEBUG
	for (int i = 0; i < NTop_Nodes; i++) {

		int nNodes = Tree[D[i].TNode.Target].DNext;

		if (D[i].TNode.Npart < 8)
			nNodes = 0;

		oprintf("DEBUG (%d:%d) Top Node %4d Target %4d nNodes %5d Npart %5d"
				" Pos %g %g %g, CoM %g %g %g \n",
				Task.Rank,Task.Thread_ID, i, D[i].TNode.Target, nNodes,
				D[i].TNode.Npart,
				D[i].TNode.Pos[0], D[i].TNode.Pos[1], D[i].TNode.Pos[2],
				D[i].TNode.CoM[0], D[i].TNode.CoM[1], D[i].TNode.CoM[2]);
	}
#endif

	Sig.Tree_Update = false;

	Profile("Build Gravity Tree");

	return ;
}

/*
 * Set the "TNode" part of the "D"omain unions. From here onwards the members 
 * are to be understood as a Top Node, not a bunch. also returns the first
 * particle in "ipart" and the "level" 
 */

static void transform_bunch_into_top_node(const int i, int *level, int *ipart)
{
	*ipart = D[i].Bunch.First_Part;
	*level = D[i].Bunch.Level;

#ifdef DEBUG
	Assert(D[i].Bunch.Npart < INT_MAX, "Npart %zu in Bunch %d > INT_MAX %d",
		   D[i].Bunch.Npart, i, INT_MAX);
#endif

	const int npart = D[i].Bunch.Npart;

	double px = P[*ipart].Pos[0] - Domain.Origin[0]; // construct from particle
	double py = P[*ipart].Pos[1] - Domain.Origin[1];
	double pz = P[*ipart].Pos[2] - Domain.Origin[2];

	double size = Domain.Size / (1ULL << *level);

	D[i].TNode.Npart = npart;
	D[i].TNode.Pos[0] = (floor(px/size) + 0.5) * size + Domain.Origin[0];
	D[i].TNode.Pos[1] = (floor(py/size) + 0.5) * size + Domain.Origin[1];
	D[i].TNode.Pos[2] = (floor(pz/size) + 0.5) * size + Domain.Origin[2];

	D[i].TNode.Target = -1;

	return ;
}

/*
 * Reserve memory in the "*Tree" structure in a thread safe way. Reallocates
 * = enlarges the *Tree memory if needed.
 */

static int reserve_tree_memory(const int i, const int nNeeded)
{
	if (nNeeded == 0)
		return 0;

	int first = 0;

	if (NNodes + nNeeded >= Max_Nodes) { // reserve more memory

		Max_Nodes = (Max_Nodes + nNeeded) * 1.05;

		Max_Nodes = imax(Max_Nodes, 1024);

		size_t nBytes = Max_Nodes * sizeof(*Tree);

		mprintf("(%d:%d) Increasing Tree Memory to %g MB, "
				"Max %d Nodes, Factor %g \n"
				, Task.Rank, Task.Thread_ID, nBytes/1024.0/1024.0, Max_Nodes,
				(double)Max_Nodes/Task.Npart_Total); fflush(stdout);

		Tree = Realloc(Tree, nBytes, "Tree");
	}

	first = NNodes;

	NNodes += nNeeded;

	return first;
}

/*
 * A subtree is build starting at node index "tree_start", from top node 
 * at index "tnode_idx". The first node is build by hand as a copy of the 
 * top node.
 * We use that the particles are in Peano-Hilbert order, i.e. that a 
 * particle will branch off as late as possible from the previous one. This 
 * means refining a node can be done via a split and reassignment of ipart and
 * ipart-1 until both are in seperate nodes. ipart+1 can then only fall
 * into the node of ipart, but not of ipart-1.
 * In the tree, DNext is the difference to the next sibling in the walk, if the
 * node is not opened. Opening a node is then node++. If DNext is negative, it 
 * points to Npart particles starting at ipart=-DNext-1, and the next node in 
 * line is node++. DNext=0 is only once per level, at the end of the branch. 
 * The tree saves only one particle per node, up to eight are combined in a 
 * node. This is achieved on the fly in an explicit cleaning step when a
 * particle opens a new branch. 
 * If the tree reaches the dynamic range of the 128bit Peano-Hilbert key, it 
 * dumps all particles in the level 42 node, making the algorithm effectively 
 * N^2  again. This happens at a depth, which corresponds to a distance of 
 * Domain.Size/2^42, hence only occurs with double precision positions. The 
 * Tree.Bitfield contains the level of the node and the Peano-Triplet of the 
 * node at that level. See Tree definition in gravity.h. 
 */

static int build_subtree(const int first_part, const int tnode_idx,
		const int top_level)
{
#ifdef DEBUG
	printf("DEBUG (%d:%d) Tree Build for top node=%d : "
		   "first part=%d npart=%d Tree build target=%d \n"
		,Task.Rank, Task.Thread_ID,  tnode_idx, first_part,
		D[tnode_idx].TNode.Npart, D[tnode_idx].TNode.Target); fflush(stdout);
#endif

	peanoKey last_key = create_first_node(first_part, tnode_idx, top_level);

	int nNodes = 1;

	int last_parent = 0; // last parent of last particle

	const int last_part = first_part + D[tnode_idx].TNode.Npart - 1;

	for (int ipart = first_part+1; ipart < last_part+1; ipart++) {

		peanoKey key = Reversed_Peano_Key(P[ipart].Pos);

		key >>= 3 * top_level;

		int node = 0;        // current node
		int lvl = top_level; // counts current level
		int parent = node;   // parent of current node

		bool ipart_starts_new_branch = true; // flag to remove leaf nodes

		while (lvl < N_PEANO_TRIPLETS) {

			if (particle_is_inside_node(key, lvl, node)) { // open node	

				if (tree[node].Npart == 1) { // refine

					tree[node].DNext = 0;

					int new_node = nNodes; // is a son of "node"

					create_node_from_particle(ipart-1, node, last_key, lvl+1,
																	new_node);
					nNodes++;

					last_key >>= 3;
				}

				add_particle_to_node(ipart, node); // add ipart to node

				ipart_starts_new_branch = ipart_starts_new_branch
													&& (node != last_parent);
				parent = node;

				node++; // decline
				lvl++;
				key >>= 3;

			} else { // skip to next node

				if (tree[node].DNext == 0 || node == nNodes - 1)
					break; // reached end of branch

				node += fmax(1, tree[node].DNext);
			}
		} // while (lvl < 42)

		if (lvl > N_PEANO_TRIPLETS-1) {

			P[ipart].Tree_Parent = parent;

			continue; // particles closer than PH resolution, dump in node
		}

		if (ipart_starts_new_branch || ipart == last_part)
			collapse_last_branch(node, last_parent, ipart, &nNodes);


		if (tree[node].DNext == 0)	// set DNext for internal node
			tree[node].DNext = nNodes - node;	// only delta

		int new_node = nNodes; // is a sibling of "node"

		create_node_from_particle(ipart, parent, key, lvl, new_node);

		nNodes++;

		last_key = key >> 3;
		last_parent = parent;

	} // for ipart

	nNodes = finalise_subtree(top_level, tnode_idx, nNodes);

#ifdef DEBUG
	printf("DEBUG (%d:%d) TNode %d tree done, nNodes %d, npart%g \n",
			Task.Rank,Task.Thread_ID, tnode_idx, nNodes,
			(double)nNodes/D[tnode_idx].TNode.Npart);
#endif

	return nNodes;
}

/*
 * The first node in the subtree needs special treatment and will later
 * be overwritten. Its properties are identical with the topnode from the
 * domain decomposition.
 */

static peanoKey create_first_node(const int first_part,
		const int tnode_idx, const int top_level)
{
	peanoKey key = Reversed_Peano_Key(P[first_part].Pos);

	key >>= 3 * top_level;

	create_node_from_particle(first_part, 0, key, top_level, 0);

	tree[0].Pos[0] = D[tnode_idx].TNode.Pos[0]; // correct position
	tree[0].Pos[1] = D[tnode_idx].TNode.Pos[1]; // because parent node 
	tree[0].Pos[2] = D[tnode_idx].TNode.Pos[2]; // did not exist

	tree[0].DUp = tnode_idx; // correct up pointer to lead to topnode

	node_set(TOP, 0);

	return key >> 3;
}

/*
 * Remove nodes from the bottom of the last branch so we end up with a
 * node with <= 8 particles. 
 */

static void collapse_last_branch(const int node, const int last_parent,
		const int ipart, int *nNodes)
{
	int n = -1;

	if (tree[node].Npart <= 8)
		n = node;
	else if (tree[last_parent].Npart <= 8)
		n = last_parent;
	else
		return ;


	int nZero = *nNodes - n - 1;

	tree[n].DNext = -ipart + tree[n].Npart - 1;

	*nNodes -= nZero;

	size_t nBytes = nZero * sizeof(*tree);

	memset(&tree[*nNodes], 0, nBytes);

	int first = -(tree[n].DNext + 1); // correct parent pointer
	int last = first + tree[n].Npart;

	for (int jpart = first; jpart < last; jpart++)
		P[jpart].Tree_Parent = n;

	return ;
}

/* 
 * Do some final operations on the node contents.
 * Copy the first node into the top node. If tree contains less than 8 
 * particles, return 0. 
 * After the build, some inner DNext pointers are 0. This is corrected 
 * setting these pointers and closing the P-H curve through the sub tree.
 * Remove the top node of the subtree, as it is identical with D.TNode.
 */

static int finalise_subtree(const int top_level, const int tnode_idx,
							int nNodes)
{
	for (int i = 0; i < nNodes; i++) { // final ops on all sub tree nodes

		tree[i].CoM[0] /= tree[i].Mass;
		tree[i].CoM[1] /= tree[i].Mass;
		tree[i].CoM[2] /= tree[i].Mass;
	}

	D[tnode_idx].TNode.Mass = tree[0].Mass; // copy first node to top node
	D[tnode_idx].TNode.CoM[0] = tree[0].CoM[0];
	D[tnode_idx].TNode.CoM[1] = tree[0].CoM[1];
	D[tnode_idx].TNode.CoM[2] = tree[0].CoM[2];

	if (tree[0].Npart <= 8) { // too small, save only topnode

		memset(tree, 0, nNodes * sizeof(*tree));

		return 0;
	}

	node_set(TOP, nNodes); // add a zero node at the end to terminate tree walk
	tree[nNodes].Mass = 1;

	nNodes++;

	tree[0].DNext = nNodes-1;  // correct internal DNext

	int stack[N_PEANO_TRIPLETS + 1] = { 0 };
	int lowest = top_level;

	for (int i = 1; i < nNodes; i++) {

		int lvl = tree[i].Bitfield & 0x3FUL; // bit 0-5

		while (lvl <= lowest) {

			int node = stack[lowest];

			if (node > 0)
				tree[node].DNext = i - node;

			stack[lowest] = 0;

			lowest--;
		}

		if (tree[i].DNext == 0) {

			stack[lvl] = i;

			lowest = lvl;
		}

	} // for i

	return nNodes;
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

	return (node_triplet == part_triplet);
}

/*
 * We always add nodes at the end of the subtree. Particles are named 
 * negative and offset by one to leave DNext=0 indicating unset. We assume
 * the peano key has reversed triplet order and the least significant 3 bits 
 * carry the triplet at level "lvl".
 */

static inline void create_node_from_particle(const int ipart,const int parent,
											 const peanoKey key, const int lvl,
											 const int node)
{
	tree[node].DNext = -ipart - 1;

	int keyfragment = (key & 0x7) << 6;

	tree[node].Bitfield = lvl | keyfragment | (1UL << 9);

	const int sign[3] = { -1 + 2 * (P[ipart].Pos[0] > tree[parent].Pos[0]),
					      -1 + 2 * (P[ipart].Pos[1] > tree[parent].Pos[1]),
		                  -1 + 2 * (P[ipart].Pos[2] > tree[parent].Pos[2]) };

	Float size = Domain.Size / (1 << lvl);

	tree[node].Pos[0] = tree[parent].Pos[0] + sign[0] * size * 0.5;
	tree[node].Pos[1] = tree[parent].Pos[1] + sign[1] * size * 0.5;
	tree[node].Pos[2] = tree[parent].Pos[2] + sign[2] * size * 0.5;

	tree[node].DUp = node - parent;

	P[ipart].Tree_Parent = node;

	add_particle_to_node(ipart, node);

	return ;
}

static inline void add_particle_to_node(const int ipart, const int node)
{
	tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;

	tree[node].Mass += P[ipart].Mass;

	tree[node].Npart++;

	return ;
}

/*
 * These functions set/extract bits from the bitmask to control
 * the tree construction and walk. They will be inlined by the compiler and
 * cost only a few cycles.  
 */

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

/*
 * this tests a sub tree for consistency, by explicitely walking the whole
 * branch of every node, collecting all contained particles and showing the
 * results alongside the values stored in the node.
 */

void test_gravity_tree(const int nNodes)
{
	for (int node = 0; node < nNodes; node++) {

		int lvl = Level(node);

		int n = node + 1;

		float mass = 0;
		float com[3] = { 0, 0, 0};
		int npart = 0;
		int nout = 0;

		double nSize = Domain.Size / (float)(1ULL << lvl);

		while (Level(n) > lvl ) { // internal node

			if (tree[n].DNext < 0) {

				int first = -tree[n].DNext - 1;
				int last = first + tree[n].Npart;

				for (int jpart = first; jpart < last; jpart++ ) {

					npart++;

					mass += P[jpart].Mass;

					com[0] += P[jpart].Pos[0] * P[jpart].Mass;
					com[1] += P[jpart].Pos[1] * P[jpart].Mass;
					com[2] += P[jpart].Pos[2] * P[jpart].Mass;

					float dx = fabs(P[jpart].Pos[0] - tree[n].Pos[0]);
					float dy = fabs(P[jpart].Pos[1] - tree[n].Pos[1]);
					float dz = fabs(P[jpart].Pos[2] - tree[n].Pos[2]);

					if (dx > nSize * 0.5)
						if (dy > nSize * 0.5)
							if (dz > nSize * 0.5)
								nout++;
				}
			}

			n++;
		}

		if (tree[node].DNext < 0) { // bundle

			int first = -tree[node].DNext - 1;
			int last = first + tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				npart++;

				mass += P[jpart].Mass;

				com[0] += P[jpart].Pos[0] * P[jpart].Mass;
				com[1] += P[jpart].Pos[1] * P[jpart].Mass;
				com[2] += P[jpart].Pos[2] * P[jpart].Mass;

				double dx = fabs(P[jpart].Pos[0] - tree[node].Pos[0]);
				double dy = fabs(P[jpart].Pos[1] - tree[node].Pos[1]);
				double dz = fabs(P[jpart].Pos[2] - tree[node].Pos[2]);

				if (dx > nSize * 0.5)
					if (dy > nSize * 0.5)
						if (dz > nSize * 0.5)
							nout++;
			}
		}

		com[0] /= mass; com[1] /= mass; com[2] /= mass;

		if (Tree[node].DNext == 0)
			printf("node %4d | m %6g =? %6g | Npart %6d =? %6d | lvl %2d "
				"| dnext %4d "
				"dup %4d nsize %8g | nout %3d | CoM %g %g %g =? %g %g %g \n",
				node, mass,tree[node].Mass, npart, tree[node].Npart, lvl,
				tree[node].DNext, tree[node].DUp, nSize, nout,
				com[0], com[1], com[2], tree[node].CoM[0],
				tree[node].CoM[1],tree[node].CoM[2] );
		Print_Int_Bits32(tree[node].Bitfield);
	}

	return ;
}


#endif // GRAVITY_TREE
