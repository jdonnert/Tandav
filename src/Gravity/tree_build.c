#include "../globals.h"
#include "../domain.h"
#include "../peano.h"
#include "../timestep.h"
#include "gravity.h"
#include "gravity_tree.h"

#ifdef GRAVITY_TREE

#define NODES_PER_PARTICLE 0.6
#define TREE_ENLARGEMENT_FACTOR 1.2

static void prepare_tree();
static bool tree_memory_is_full(const int );
static int reserve_tree_memory(const int);
static void set_tree_parent_pointers (const int);
static int build_subtree(const int, const int, const int);
static int finalise_subtree(const int, const int, int );
static inline bool particle_is_inside_node(const peanoKey,const int,const int);
static inline void add_particle_to_node(const int, const int);
static peanoKey create_first_node(const int, const int, const int);
static void collapse_last_branch(const int, const int, const int, int*);
static inline int key_fragment(const int);
static inline void node_set(const enum Tree_Bitfield, const int);
static void print_top_nodes();
static inline void create_node_from_particle(const int, const int,
											 const peanoKey, const int,
											 const int);
size_t buf_threshold = 0;
uint32_t NNodes = 0;
static int Max_Nodes = 0;
struct Tree_Node  * restrict Tree = NULL; // global pointer to all nodes
static omp_lock_t Tree_Lock; // lock global *Tree, NNodes, Max_Nodes

static struct Tree_Node * restrict tree = NULL; //  build in *Tree or *Buffer
#pragma omp threadprivate(tree)

/*
 * This builds the tree in parallel, particles are assumed PH ordered. 
 * We build the tree corresponding to a top node either in the openmp buffer, 
 * or directly inside the *Tree memory. Every access to  NNodes, 
 * Max_Nodes has to be protected by the openmp Tree_Lock.
 * A subtree is build starting from the top node in target pointer "*tree". 
 * Every subtree allocated a fixed amount of memory, if it is not build in the
 * buffer.
 * If the number of particles in that top node is <= 8, the subtree is 
 * discarded and the topnode target points directly to the particles and 
 * the tree walk will use particles directly from the topnode. 
 * At last the top nodes are broadcasted.
 */

void Gravity_Tree_Build()
{
	Profile("Build Gravity Tree");

	#pragma omp single
	NNodes = 0;

	for (;;) {

		#pragma omp single
		prepare_tree();

		#pragma omp for schedule(static,1)
		for (int i = 0; i < NTop_Nodes; i++) {

	//		if (D[i].TNode.Target < 0)  // not local
	//			continue;
	
			bool build_in_buffer = D[i].TNode.Npart < buf_threshold;

			int nReserved = ceil(D[i].TNode.Npart * NODES_PER_PARTICLE);

			if (build_in_buffer) {

				tree = Get_Thread_Safe_Buffer(Task.Buffer_Size);

			} else { // build in *Tree directly

				D[i].TNode.Target = reserve_tree_memory(nReserved);
	
				tree = &Tree[D[i].TNode.Target];
			}

			if (tree_memory_is_full(nReserved))
				continue;
			
			int nNeeded = build_subtree(D[i].TNode.First_Part, i, 
					D[i].TNode.Level);

			if (build_in_buffer) { // copy buffer to Tree

				if (tree_memory_is_full(nNeeded))
					continue;

				D[i].TNode.Target = reserve_tree_memory(nNeeded);

				size_t nBytes = nNeeded * sizeof(*Tree);

				memcpy(&Tree[D[i].TNode.Target], tree, nBytes);
			}	

			set_tree_parent_pointers(i);
			
		} // for i

		#pragma omp barrier

		if (! tree_memory_is_full(0))
			break;

		#pragma omp barrier

	} // forever

	rprintf("Tree build: %d of %d Nodes used (%g MB)\n",
			NNodes, Max_Nodes, Max_Nodes*sizeof(*Tree)/1024.0/1024);

	print_top_nodes(); // DEBUG only

	Sig.Tree_Update = false;

	Profile("Build Gravity Tree");

	return ;
}

void Setup_Gravity_Tree()
{
	omp_init_lock(&Tree_Lock); // we don't destroy this one ...

	Max_Nodes = 0.3 * Task.Npart_Total;
			
	Tree = Malloc(Max_Nodes * sizeof(*Tree), "Tree");

	buf_threshold = Task.Buffer_Size/sizeof(*Tree);

	return ;
}

static void prepare_tree()
{
	if (NNodes != 0) { // tree build aborted, increase mem

		Max_Nodes = TREE_ENLARGEMENT_FACTOR * Max_Nodes;
		
		printf("(%d:%d) Increased Tree Memory to %6.1f MB, "
			"Max %10d Nodes, ratio %4g \n", Task.Rank, Task.Thread_ID, 
			Max_Nodes * sizeof(*Tree)/1024.0/1024.0, Max_Nodes, 
			(double) Max_Nodes/Task.Npart_Total); 
	}

	Tree = Realloc(Tree, Max_Nodes * sizeof(*Tree), "Tree");
		
	memset(Tree, 0, Max_Nodes * sizeof(*Tree));

	NNodes = 0;

	return ;
}

/*
 * Reserve memory in the "*Tree" structure by increasing NNodes. We lock 
 * NNodes with Tree_Lock, so threads don't mess up the pointers. 
 */

static int reserve_tree_memory(const int nNeeded)
{
	omp_set_lock(&Tree_Lock);

	int first = NNodes;

	NNodes += nNeeded;

	omp_unset_lock(&Tree_Lock);

	return first;
}

/*
 * Check if we have enough memory left to build/copy-in the tree. If not
 * redo the whole tree build
 */

static bool tree_memory_is_full(const int nReserved)
{
	bool result = false;
	
	omp_set_lock(&Tree_Lock);

	if (NNodes + nReserved > Max_Nodes) {

		NNodes = Max_Nodes + 1; // trigger another tree build
	
		result = true;
	}

	omp_unset_lock(&Tree_Lock);

	return result;
}

/* 
 * Correct the Tree_Parent pointers in P in case we build in the buffer
 * which always starts at 0. If the top node doesn't contain a tree, make a 
 * pointer to the top node instead.
 */

static void set_tree_parent_pointers (const int i)
{
	const int first_part = D[i].TNode.First_Part;
	const int last_part = first_part + D[i].TNode.Npart;

	if (D[i].TNode.Target > 0) { // correct particle parent pointer
	
		for (int ipart = first_part; ipart < last_part; ipart++)
			P[ipart].Tree_Parent += D[i].TNode.Target;
		
	} else if (D[i].TNode.Target < 0) {
	
		for (int ipart = first_part; ipart < last_part; ipart++)
			P[ipart].Tree_Parent = -i - 1; // top node w/o tree
	}
	return ;
}

/*
 * A subtree is build starting in *tree, from top node at index "tnode_idx". 
 * The first node is build by hand from the existing top node.
 * We use that the particles are in Peano-Hilbert order, i.e. that a 
 * particle will branch off as late as possible from the previous one. This 
 * means refining a node can be done via a split and reassignment of ipart and
 * ipart-1 until both are in seperate nodes. ipart+1 can then only fall
 * into the node of ipart, but not of ipart-1.
 * In the tree, DNext is the difference to the next sibling in the walk, if the
 * node is not opened. Opening a node is then node++. If DNext is negative, it 
 * points to Npart particles starting at ipart=-DNext-1, and the next node in 
 * line is node++. DNext=0 is only once, at the end of the branch. 
 * The tree saves only one particle per node, up to eight are combined in a 
 * node. This is achieved on the fly in an explicit cleaning step when a
 * particle opens a new branch. 
 * If the tree reaches the dynamic range of the 128bit Peano-Hilbert key, it 
 * dumps all particles in the level 42 node, making the algorithm effectively 
 * N^2  again. This happens at a depth, which corresponds to a distance of 
 * Domain.Size/2^42, hence only occurs with double precision positions. The 
 * Tree.Bitfield contains the level of the node and the Peano-Triplet of the 
 * node at that level. See *Tree definition in gravity.h. 
 */

static int build_subtree(const int first_part, const int tnode_idx,
		const int top_level)
{
#ifdef DEBUG
	printf("DEBUG (%d:%d) Tree Build for top node=%d : "
		   "first part=%d npart=%d Tree build target=%d  \n"
		,Task.Rank, Task.Thread_ID,  tnode_idx, first_part,
		D[tnode_idx].TNode.Npart, D[tnode_idx].TNode.Target);
	fflush(stdout);
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
	printf("DEBUG (%d:%d) TNode %d tree done, nNodes %d, npart/nNodes %g \n",
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
 * Copy the first node into the top node. If the subtree contains less than 8 
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

	if (tree[0].Npart <= 8) { // too small, save only topnode, return empty

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

	return node_triplet == part_triplet;
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

static void print_top_nodes()
{
#ifdef DEBUG
	#pragma omp single
	for (int i = 0; i < NTop_Nodes; i++) 
		printf("%d Target=%d Level=%d Npart=%d Pos=%g %g %g, "
				"Mass=%g CoM=%g %g %g Dp=%g %g %g \n",i, 
				D[i].TNode.Target,D[i].TNode.Level, D[i].TNode.Npart,
				D[i].TNode.Pos[0],D[i].TNode.Pos[1],D[i].TNode.Pos[2],
				D[i].TNode.Mass,
				D[i].TNode.CoM[0],D[i].TNode.CoM[1],D[i].TNode.CoM[2],
				D[i].TNode.Dp[0],D[i].TNode.Dp[1],D[i].TNode.Dp[2]);

#endif
	
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
