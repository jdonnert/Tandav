#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#ifdef GRAVITY_TREE

#define NODES_PER_PARTICLE 0.7 

struct Tree_Node *Tree = NULL;
size_t NNodes = 0;
size_t Max_Nodes = 0;

static inline int level(const int node); // bitmask functions
static inline int key_fragment(const int node);

static int build_subtree(const int, const int, const int, const int);
static inline void add_node(const int, const int, const peanoKey, 
		const int, int *);
static inline bool particle_is_inside_node(const peanoKey *, const int node);
static inline void add_particle_to_node(const int ipart, const int node);

void gravity_tree_init();

static void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--) {
		printf("%llu", (val & (1ULL << i) ) >> i);
		if (i % 3 == 0 && i != 0)
			printf(".");
	}
	printf("\n");fflush(stdout);

	return ;
}

/*
 * This builds the tree. 
 * First the tree in every local node is build starting from the bunchnode. 
 * In the global tree memory, all bunchnodes are seperated by npart nodes.
 * Then the tree from the top nodes is constructed and the local nodes 
 * are moved in blocks.
 *
 */

void Gravity_Tree_Build()
{
	Profile("Build Gravity Tree");

	gravity_tree_init();

	#pragma omp single
	{
	// build_top_tree();
	// compact_TREE();
	// do_final_operations();
	NNodes = build_subtree(0, Task.Npart_Total, 0, 0);

	printf("\nTree build: zu nodes for %d particles\n", 
			NNodes, Task.Npart_Total );

	} // omp single
	
	Profile("Build Gravity Tree");


	return ;
}

static void build_top_tree()
{

	return ;
}

static void compact_tree()
{

	return ;
}

static void finalise_tree()
{
	#pragma omp for
	for (int i = 0; i < NNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;

		Tree[i].Vel_CoM[0] /= Tree[i].Mass;
		Tree[i].Vel_CoM[1] /= Tree[i].Mass;
		Tree[i].Vel_CoM[2] /= Tree[i].Mass;
	}
}


/*
 * A subtree is build starting at node index top.
 * In the tree, DNext points to the difference to next index in the walk, 
 * if the node is not opened. Opening a node is then node++. If DNext is 
 * negative it points to Npart particles starting at ipart=-DNext-1, and the 
 * next node in line is node++. DNext=0 is only at the root node.
 * The tree saves only one particle per node, up to eight are combined in
 * the node. This is done on the fly in an explicit cleaning step, when a 
 * particle opens a new branch. 
 * The Bitfield contains the level of the node and the Peano-Triplet of the
 * node at that level. See Tree definition.
 * After the initial build, some DNext pointers are 0 and a non-recursive tree
 * walk is performed to set these pointers and close the P-H curve through the
 * tree.
 */

static int build_subtree(const int istart, const int npart, const int offset, 
		const int top)
{
#ifdef DEBUG
	printf("DEBUG: (%d:%d) Sub-Tree Build istart=%d npart=%d offs=%d top=%d \n"
			, Task.Rank, Task.Thread_ID, istart, npart, offset, top);
#endif

	int nNodes = 0; // local in this subtree

	Float px = (P[0].Pos[0] - Domain.Origin[0]) / Domain.Size;
	Float py = (P[0].Pos[1] - Domain.Origin[1]) / Domain.Size;
	Float pz = (P[0].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
	peanoKey pkey[2] = { 0 };
	pkey[0] = Peano_Key(px, py, pz, &pkey[1]);
	
	__uint128_t oldkey = 0;
	memcpy(*pkey[0], &oldkey, 32);
	add_node(istart, offset, pkey[0], top, &nNodes); // add the first node by hand

	Tree[0].Pos[0] = Domain.Center[0];
	Tree[0].Pos[1] = Domain.Center[1];
	Tree[0].Pos[2] = Domain.Center[2];

	int last_parent = offset;		// last parent of last particle

	for (int ipart = istart+1; ipart < Task.Npart_Total; ipart++) {

		Float px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		Float py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		Float pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		peanoKey pkey[2] = { 0 };
		pkey[0] = Peano_Key(px, py, pz, &pkey[1]);

		int node = offset;  		// current node
		int lvl = top;				// counts current level
		int parent = node;			// parent of current node

		bool ipart_starts_new_branch = true; // flag to remove leaf nodes
		
		__uint128_t key = 0;
		memcpy(*pkey[0], &key, 32);
		
		for (;;) {

			if (particle_is_inside_node(key, node)) { // open node
			
				if (Tree[node].DNext == -ipart) { // points to ipart-1, refine 
				
					Tree[node].DNext = 0; // mark DNext free
					
					if (lvl+1 < 22) // add daughter
						add_node(ipart-1,node, pkey[0], lvl+1,&nNodes);  
					else
						add_node(ipart-1,node, pkey[1], lvl+1,&nNodes)
				}  
				
				add_particle_to_node(ipart, node); // add ipart to node

				ipart_starts_new_branch &= (node != last_parent);

				parent = node;

				node++; // decline into node

				lvl++;

			} else { // skip node
				
				if (Tree[node].DNext == 0 || node == nNodes - 1)   
					break; // reached end of my branch or whole tree
				
				node += fmax(1, Tree[node].DNext);
			}

			key <<= 3;

		} // for (;;)

		if (ipart_starts_new_branch) { // now remove leaf particle nodes

			int n = 0;

			if (Tree[node].Npart <= 8) 
				n = node;			
			else if (Tree[last_parent].Npart <= 8) 
				n = last_parent;
			
			if (n != 0) {

				Tree[n].DNext = -ipart + Tree[n].Npart - 1;
					
				int nZero = nNodes - n - 1;

				memset(&Tree[n+1], 0, nZero*sizeof(*Tree));
				
				nNodes = n + 1;
			}	
		}
	
		if (Tree[node].DNext == 0) // set internal next node
			Tree[node].DNext = nNodes - node; // only delta

		if (lvl < 22) // add sibling to current node
			add_node(ipart, parent, pkey[0], lvl, &nNodes); 
		else
			add_node(ipart, parent, pkey[1], lvl, &nNodes); 
	
		last_parent = parent;
		oldkey = key;
	
	} // for ipart

	//#pragma omp single nowait // correct for last DNext to point upwards
	{

	Tree[top].DNext = 0; 

	int stack[63] = { 0 }; 
	int lowest = 0;

	for (int i = 1; i < nNodes; i++) {
		
		int lvl = level(i);

		while (lvl <= lowest) {

			int node = stack[lowest];
	
			if (node > 0)
				Tree[node].DNext = i - node;

			stack[lowest] = 0;

			lowest--;
		} 
		
		if (Tree[i].DNext == 0) { 
			
			stack[lvl] = i;
			
			lowest = lvl;
		}
		
	} // for
	
	} // omp single nowait

	for (int i = 0; i < nNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;

		Tree[i].Vel_CoM[0] /= Tree[i].Mass;
		Tree[i].Vel_CoM[1] /= Tree[i].Mass;
		Tree[i].Vel_CoM[2] /= Tree[i].Mass;
	}

	return nNodes;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal. If the tree is deeper than the number of 
 * bit triplets in the peano key (i.e. level>21), we always return 
 * "false", causing the particles to be added as leaves at the end of the 
 * tree in random order. The tree walk then has to look at all of them
 * which make this N^2.
 */

static inline bool particle_is_inside_node(const peanoKey pkey,
		const int node)
{
	const int lvl = level(node);

	const int shift = 63 - 3 * lvl; // bits 0-5

	const uint64_t part_mask = 0x07ULL << shift; // shift first 3 bits

	uint64_t part_triplet = (key & part_mask) >> shift;

	uint64_t node_triplet = key_fragment(node); 

	return (node_triplet == part_triplet); 
}

/*
 * We always add nodes at the end of the tree. Particles are negative and 
 * offset by one to leave DNext=0 indicating unset.
 */

static inline void add_node(const int ipart, const int parent, 
		const peanoKey key, const int lvl, int *nNodes)
{
	const int node = (*nNodes)++;

	Assert(node < Max_Nodes, "Too many nodes (%d), increase "
			"NODES_PER_PARTICLE=%d", nNodes, NODES_PER_PARTICLE);

	Tree[node].DNext = -ipart - 1;

	int shift = 63 - 3*(lvl % 22);

	uint32_t keyfragment = ((key >> shift) & 0x7) << 6; 

	Tree[node].Bitfield = lvl | keyfragment;
	
	const int sign[3] = { -1 + 2 * (P[ipart].Pos[0] > Tree[parent].Pos[0]),
	 			     	  -1 + 2 * (P[ipart].Pos[1] > Tree[parent].Pos[1]),
	 			          -1 + 2 * (P[ipart].Pos[2] > Tree[parent].Pos[2])}; 
	
	Float size = Domain.Size / (1 << lvl);

	Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size * 0.5;
	Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size * 0.5;
	Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size * 0.5;

	//Tree[node].DUp = node - parent;

	add_particle_to_node(ipart, node); 

	return ;
}

static inline void add_particle_to_node(const int ipart, const int node)
{
	Tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	Tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	Tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;

	Tree[node].Vel_CoM[0] += P[ipart].Vel[0] * P[ipart].Mass;
	Tree[node].Vel_CoM[1] += P[ipart].Vel[1] * P[ipart].Mass;
	Tree[node].Vel_CoM[2] += P[ipart].Vel[2] * P[ipart].Mass;
	
	Tree[node].Mass += P[ipart].Mass;

	Tree[node].Npart++;
	
	return ;
}

static inline void add_bunch_to_node(const int i, const int node)
{
	Tree[node].CoM[0] += Bunch[i].CoM[0];
	Tree[node].CoM[1] += Bunch[i].CoM[1];
	Tree[node].CoM[2] += Bunch[i].CoM[2];

	Tree[node].Mass += Bunch[i].Mass;

	Tree[node].Npart += Bunch[i].Npart;

	return ;
}

/*
 * These functions set/extract bits from the bitmask to control
 * the tree construction and walk. They will be inlined by the compiler and
 * cost only a few cycles. We store bitmasks separately to ensure the
 * correct type of the constant.
 */

static inline int level(const int node)
{
	const uint32_t bitmask = 0x3F;

	return Tree[node].Bitfield & bitmask; // return but 0-5
}

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (Tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

void gravity_tree_init()
{
	#pragma omp single
	{

	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	size_t nBytes = Max_Nodes * sizeof(*Tree);

	if (Tree == NULL)
		Tree = Malloc(nBytes, "Tree");
	
	memset(Tree, 0, nBytes);

	NNodes = 0;

	Tree[0].Pos[0] = Domain.Center[0];
	Tree[0].Pos[1] = Domain.Center[1];
	Tree[0].Pos[2] = Domain.Center[2];

	} // omp single

	return ;
}

void Gravity_Tree_Update()
{

	return ;
}

void test_gravity_tree(const int nNodes)
{
	for (int node = 0; node < nNodes; node++) {
	
		int lvl = level(node);
	
		int n = node + 1;

		float mass = 0;
		int npart = 0;
		int nout = 0;

		double nSize = Domain.Size / (float)(1ULL << lvl);

		if (Tree[node].DNext < 0)
			continue;

		while (level(n) > lvl) {

			if (Tree[n].DNext < 0) {
			
				int first = -Tree[n].DNext - 1;
				int last = first + Tree[n].Npart;

				for (int jpart = first; jpart < last; jpart++ ) {

				 	npart++;

					mass += P[jpart].Mass;
	
					float dx = fabs(P[jpart].Pos[0] - Tree[n].Pos[0]);
					float dy = fabs(P[jpart].Pos[1] - Tree[n].Pos[1]);
					float dz = fabs(P[jpart].Pos[2] - Tree[n].Pos[2]);

					if (dx > nSize * 0.5) 
						if (dy > nSize * 0.5) 
							if (dz > nSize * 0.5)
								nout++;
				}
			}

			n++; 
		}

		printf("%d m=%g,%g N=%d,%d nsize=%g nout=%d\n", node, mass, 
				Tree[node].Mass, npart, Tree[node].Npart, 	nSize, nout);
	}

	return ;
}


#endif // GRAVITY_TREE
