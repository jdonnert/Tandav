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
static inline void add_node(const int, const int, const int, int *);
static inline bool particle_is_inside_node(const int ipart, const int node);
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

void Build_Gravity_Tree()
{
	Profile("Build Gravity Tree");

	gravity_tree_init();

	#pragma omp single
	{
	
	NNodes = build_subtree(0, Task.Npart_Total, 0, 0);

	printf("Finished tree build, %d nodes for %d particles \n", 
			NNodes, Task.Npart_Total);
	}
	Profile("Build Gravity Tree");


	return ;
}

/*
 * A subtree is build starting at node index top.
 * In the tree, DNext points to the difference to next index in the walk, 
 * if the node is not opened. Opening a node is then node++. If DNext is 
 * negative it points to Nprat particles starting at ipart=-DNext-1, and the 
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

	int nNodes = 0;

	add_node(istart, offset, top, &nNodes); 	// add the first node by hand

	Tree[0].Pos[0] = Tree[0].Pos[1] = Tree[0].Pos[2] = 0;

	int last_parent = offset;		// last parent of last particle

	for (int ipart = istart+1; ipart < Task.Npart_Total; ipart++) {
		
		int node = offset;  		// current node
		int lvl = top;				// counts current level
		int parent = node;			// parent of current node

		bool ipart_starts_new_branch = true; // flag to remove leaf nodes

		for (;;) {

			if (particle_is_inside_node(ipart, node)) { // open node
			
				if (Tree[node].DNext == -ipart) { // points to ipart-1, refine 
				
					Tree[node].DNext = 0; // mark DNext free
					
					add_node(ipart-1, node, lvl+1, &nNodes); // add daughter 
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

		add_node(ipart, parent, lvl, &nNodes); // add sibling to current node
	
		last_parent = parent;
	
	} // for ipart

	#pragma omp single nowait // correct for last DNext to point upwards
	{

	Tree[top].DNext = 0;

	int stack[63] = { 0 }; // can't go deeper than 63
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

	#pragma omp for
	for (int i = 0; i < nNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;
	}
/*
for (int n = 0; n < nNodes; n++) {
	//if (level(n) == 2) {
printf("TEST n=%d np=%d next=%d  mass=%g level=%d  \n", 
n,  Tree[n].Npart, Tree[n].DNext, Tree[n].Mass, level(n));
//print_int_bits64(Tree[n].Bitfield);
	//}
}
printf("\n");

for (int ipart = 0; ipart < 56; ipart++) { 
printf("%d ", ipart); print_int_bits64(P[ipart].Peanokey);}
exit(0); */

//for (int z = 0; z < nNodes; z++) printf("TREE n=%d np=%d next=%d  mass=%g level=%d  \n", z,  Tree[z].Npart, Tree[z].DNext, Tree[z].Mass, level(z));

	for (int node = 0; node < nNodes; node++) {
	
		int lvl = level(node);
	
		int n = node + 1;

		double mass = 0;
		int npart = 0;
		int nout = 0;

		float nSize = Domain.Size / (float)(1ULL << lvl);

		while (level(n) > lvl) {
//printf("n=%d lvl=%d np=%d \n", n, level(n), Tree[n].Npart);


			if (Tree[n].DNext < 0) {
			
				int first = -Tree[n].DNext - 1;
				int last = first + Tree[n].Npart;

				for (int jpart = first; jpart < last; jpart++ ) {

				 	npart++;

					mass += P[jpart].Mass;
//printf("    jpart=%d first=%d last=%d np=%d lvl=%d \n", jpart, first, last, Tree[n].Npart, level(n));
	
					float dx = fabs(P[jpart].Pos[0] - Tree[n].Pos[0]);
					float dy = fabs(P[jpart].Pos[1] - Tree[n].Pos[1]);
					float dz = fabs(P[jpart].Pos[2] - Tree[n].Pos[2]);

					if (dx > nSize * 0.5)
					if (dy > nSize * 0.5)
					if (dz > nSize * 0.5)
						nout++;
				}
			}

			n++; //= fmax(1, Tree[n].DNext);
		}

printf("%d m=%g,%g N=%d,%d nsize=%g nout=%d\n", 
node, mass, Tree[node].Mass, npart, Tree[node].Npart, 	nSize, nout);
	
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

static inline bool particle_is_inside_node(const int ipart, const int node)
{
	const int lvl = level(node);

	const int shift = 63 - 3 * lvl; // bits 0-5

	const uint64_t part_mask = 0x07ULL << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].Peanokey & part_mask) >> shift;

	uint64_t node_triplet = key_fragment(node); 

	return (node_triplet == part_triplet) && (lvl < 22); // branch free
}

/*
 * We always add nodes at the end of the tree. Particles are negative and 
 * offset by one to leave DNext=0 indicating unset.
 */

static inline void add_node(const int ipart, const int parent, const int lvl, 
		int *nNodes)
{
	const int node = (*nNodes)++;

	Tree[node].DNext = -ipart - 1;
	
	int shift = 63 - 3*lvl - 6 ; // leave 6 bits space for level in bitfield
	
	uint32_t keyfragment = (P[ipart].Peanokey >> shift) & 0x1C0ULL;

	Tree[node].Bitfield = lvl | keyfragment; 
	
	int sign[3] = { 1 - 2 * (int)(P[ipart].Pos[0] < Tree[parent].Pos[0]),
	 			    1 - 2 * (int)(P[ipart].Pos[1] < Tree[parent].Pos[1]),
	 			    1 - 2 * (int)(P[ipart].Pos[2] < Tree[parent].Pos[2])}; 
	
	Float size = Domain.Size / (1 << (lvl));

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

	Tree[0].Pos[0] = Tree[0].Pos[1] = Tree[0].Pos[2] = 0;

	} // omp single

	return ;
}


#endif // GRAVITY_TREE
