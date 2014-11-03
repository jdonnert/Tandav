#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#define NODES_PER_PARTICLE 1.0 

struct Tree_Node {
	uint32_t Bitfield; // bit 0-5:level, 6-8:key, 9: DNext flag, 10-31:free
	int DNext;		   // this is the DELTA to the next node; or -ipart
	Float CoM[3];
	float Mass;
	//int DUp;
	int Npart;
} *Tree = NULL;

static size_t NNodes = 0;
static size_t Max_Nodes = 0;

static inline int level(const int node); // bitmask functions
static inline int key_fragment(const int node);
static inline Float node_size(const int node);
static inline void add_particle_to_bitfield(const int node);

static inline void add_node(const int ipart, const int node, const int lvl);
static inline bool particle_is_inside_node(const int ipart, const int node);
static inline void add_particle_to_node(const int ipart, const int node);

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
 *
 */

void Build_Tree()
{
	Profile("Build Gravity Tree");

	Init_Tree();

	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/
 
 	/* Here the top tree is complete, now build local nodes downwards */

	add_node(0, 0, 0); // add the first particle by hand
	Tree[0].Bitfield &= ~0x1FF; // clear first 9 bits

	int last_parent = 0; 		// parent of last particle

	for (int ipart = 1; ipart < Task.Npart_Total; ipart++) {

		int node = 0;  			// current node
		int lvl = 0;			// counts current level
		int parent = 0;			// parent of current node
		bool is_new_branch = true;	// flag to delete single particle nodes

		for (;;) {

if (level(node) != lvl) {
printf("ERROR LEVEL! %d: %d!=%d \n",node, level(node), lvl); goto out;}

			if (particle_is_inside_node(ipart, node)) { // open node
			
				if (Tree[node].DNext == -ipart) { // refine 
				
					int jpart = ipart-1; // node has to contain last particle 
					
					Tree[node].DNext = 0; // mark DNext free

					add_node(jpart, node, lvl+1); // add daughter to node 
				}  
				
				add_particle_to_node(ipart, node); // ipart to current node

				if (node == last_parent)
					is_new_branch = false;

				parent = node;
				node++; // decline into node containing jpart
				lvl++;

			} else { // skip

				if (Tree[node].DNext == 0 || node == NNodes - 1)   
					break; // reached end of my branch or whole tree
				
				if (Tree[node].DNext < 0)  // internal leaf
					node++; // might need to split its sibling 
				else 
					node += Tree[node].DNext; // skip whole subtree
			}

		} // for (;;)

		if (is_new_branch) { // remove leaf particle nodes from tree

			int n = NNodes-1; 	// node counter
			int np = 0;			// particle counter
			
			while (Tree[n].DNext < 0) { // walk backwards
			
				n--;
				np++;
			}

			if (np == 0)
				break;

			Tree[n].DNext = Tree[n+1].DNext;
			
			NNodes = n + 1; // remove leaf nodes
			memset(&Tree[n+1], 0, np*sizeof(*Tree));
		}
	
		if (Tree[node].DNext == 0) // set internal next node
			Tree[node].DNext = NNodes - node; // only delta

		add_node(ipart, parent, lvl); // add a sibling to current node
	
		last_parent = parent;

	} // for ipart
out:;

	#pragma omp single nowait // correct DNext=0 to point upwards
	{

	Tree[0].DNext = -1;

	int stack[25] = { 0 };
	int lowest = 0;

	for (int i = 1; i < NNodes; i++) {
		
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
	
	} // omp single

for (int n=0; n<Task.Npart_Total; n++) {
	if (level(n) == 2) {
printf("TEST n=%d np=%d next=%d  mass=%g level=%d  \n", 
n,  Tree[n].Npart, Tree[n].DNext, Tree[n].Mass, level(n));
print_int_bits64(Tree[n].Bitfield);
	}
}
printf("\n");

for (int ipart = 0; ipart < 10; ipart++) { 
printf("%d ", ipart); print_int_bits64(P[ipart].Peanokey);}

	#pragma omp for
	for (int i = 0; i < NNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;
	}

	rprintf("Finished tree build, %d nodes for %d particles", 
			NNodes, Task.Npart_Total);

	Profile("Build Gravity Tree");
exit(0);
	return ;
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

static inline void add_node(const int ipart, const int parent, const int lvl)
{
	const int node = NNodes++;

	Tree[node].DNext = -ipart - 1;
	
	int shift = 63 - 3*lvl - 6 ; // leave 6 bits space for level in bitfield
	
	uint32_t keyfragment = (P[ipart].Peanokey >> shift) & 0x1C0ULL;

	Tree[node].Bitfield = lvl | keyfragment; 

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

static inline Float node_size(const int node)
{
	int lvl = level(node);

	return Domain.Size / ((Float) (1UL << lvl)); // Domain.Size/2^level
}


void Init_Tree()
{
	#pragma omp single
	{

	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	size_t nBytes = Max_Nodes * sizeof(*Tree);

	if (Tree == NULL)
		Tree = Malloc(nBytes, "Tree");
	
	memset(Tree, 0, nBytes);

	Tree[0].CoM[0] = 0; // will hold global CoM of sim
	Tree[0].CoM[1] = 0;
	Tree[0].CoM[2] = 0;

	Tree[0].Npart = 0;
	//Tree[0].DUp = -1;
	Tree[0].Mass = 0;
	Tree[0].DNext = -1;
	Tree[0].Bitfield = 0;

	NNodes = 0;

	} // omp single

	return ;
}

/*
 * This function computes the gravitational acceleration by walking the tree  
 */

void Accel_Gravity_Tree(const int ipart, double *acc, double *pot)
{

	if (Sig.Fullstep) // need to build new one
		Free(Tree);

	return ;
}
