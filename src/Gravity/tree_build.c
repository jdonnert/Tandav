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
 * if the node is not opened, opening is node++. If DNext is negative it 
 * points to particles with ipart=-DNext-1, and the next node in line 
 * is node++. 
 * The tree saves only one particle per leaf-node. Other particles in the 
 * same leaf node have to be inferred from Npart during the walk, to 
 * safe memory.
 * The Bitfield contains the level of the node and the Peano-Triplet of the
 * node at that level. 
 * After the initial build, some DNext pointers are 0 and a non-recursive tree
 * walk is performed to set these pointers and close the P-H curve through the
 * tree.
 */

static int build_subtree(const int start, const int npart, const int offset, 
		const int top)
{
	int nNodes = 0;

	add_node(start, offset, top, &nNodes); 	// add the first particle by hand

printf("Tree Build Local: istart=%d npart=%d offs=%d top=%d nNodes=%d\n", 
start, npart, offset, top, nNodes);

	int last_parent = offset;		// parent of last particle

	for (int ipart = start+1; ipart < npart; ipart++) {
		
		int node = offset;  		// current node
		int lvl = top;				// counts current level
		int parent = node;			// parent of current node
		bool is_new_branch = true;	// flag to delete single particle nodes

		for (;;) {

if (level(node) != lvl) {
printf("ERROR LEVEL! %d: %d!=%d \n",node, level(node), lvl); goto out;}

			if (particle_is_inside_node(ipart, node)) { // open node
			
				if (Tree[node].DNext == -ipart) { // refine 
				
					Tree[node].DNext = 0; // mark DNext free

					add_node(ipart-1, node, lvl+1, &nNodes); // add daughter 
				}  
				
				add_particle_to_node(ipart, node); // add ipart to node

				is_new_branch = is_new_branch & (node != last_parent);

				parent = node;

				node++; // decline into node containing jpart

				lvl++;

			} else { // skip

				if (Tree[node].DNext == 0 || node == nNodes - 1)   
					break; // reached end of my branch or whole tree
				
				node += fmax(1, Tree[node].DNext);
			}

		} // for (;;)

		if (is_new_branch) { // remove leaf particle nodes from tree

			int n = nNodes-1; 	// node counter
			int np = 0;			// particle counter
			
			while (Tree[n].DNext < 0) { // walk backwards
			
				n--;
				np++;
			}

			if (np == 0)
				break;

			Tree[n].DNext = Tree[n+1].DNext;
			
			nNodes = n + 1; // remove leaf nodes

			memset(&Tree[n+1], 0, np*sizeof(*Tree));
		}
	
		if (Tree[node].DNext == 0) // set internal next node
			Tree[node].DNext = nNodes - node; // only delta

		add_node(ipart, parent, lvl, &nNodes); // add sibling to current node
	
		last_parent = parent;

	} // for ipart
out:;

	#pragma omp single nowait // correct DNext=0 to point upwards
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

	#pragma omp for
	for (int i = 0; i < nNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;
	}
/*
for (int n = 0; n < Task.Npart_Total; n++) {
	if (level(n) == 2) {
printf("TEST n=%d np=%d next=%d  mass=%g level=%d  \n", 
n,  Tree[n].Npart, Tree[n].DNext, Tree[n].Mass, level(n));
print_int_bits64(Tree[n].Bitfield);
	}
}
printf("\n");

for (int ipart = 0; ipart < 10; ipart++) { 
printf("%d ", ipart); print_int_bits64(P[ipart].Peanokey);}
*/
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

	//Tree[node].DUp = node - parent;
	
	int sign[3] = { -1 * (int)(P[ipart].Pos[0] < Tree[parent].Pos[0]),
	 			    -1 * (int)(P[ipart].Pos[1] < Tree[parent].Pos[1]),
	 			    -1 * (int)(P[ipart].Pos[2] < Tree[parent].Pos[2])}; 
	
	Float size = Domain.Size / (1 << lvl+1);

	Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size * 0.5;
	Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size * 0.5;
	Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size * 0.5;

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

	Tree[0].Pos[0] = Tree[0].Pos[1] = Tree[0].Pos[2] = 1;

	} // omp single

	return ;
}


#endif // GRAVITY_TREE
