/*
 * This builds the tree. Nodes are ordered according to their Peano-Hilbert
 * key (triplet) on that level. Outside of the key range, nodes are shuffled
 * and the treewalk has to open them all. However this will be close to the
 * precision limit of floats: 2^-21.
 * We store level and key triplet in a 32 bit mask. The peano order implies 
 * that only non-empty nodes are constructed and we do not need a "down" 
 * pointer, which is just a node++ now. It also means that new nodes are 
 * always added at the end of the tree during construction and the tree walk
 * has best possible cache coherency. 
 * To make the subtree movable, the "next" & "up" pointers just give the 
 * difference between current and next/parent node.
 * For particles closer than box/2^21, i.e. with the same 63bit peano key, 
 * the nodes are just randomly inserted. This implies that the tree walk has
 * to open all of these later and becomes ~N^2 for this subtree.
 */

#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#define NODES_PER_PARTICLE 5 // Springel 2005

struct Tree_Node {
	uint32_t Bitfield; // bit 0-5:level, 6-8:key, 9-32:free
	int Npart;
	Float CoM[3];
	float Mass;
	int DNext;		 // this is the DELTA to the next node
	int DUp;
} *Tree = NULL;


static size_t NNodes = 0;
static size_t Max_Nodes = 0;

static inline int level(const int node);
static inline int key_fragment(const int node);
static inline Float node_size(const int node);

static inline void add_node(const int ipart, const int node);
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


void Build_Tree()
{
	Init_Tree();

	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/
 
 	/* Here the top tree is complete, now build local nodes downwards */

for (int n=0; n<NNodes; n++)
	printf("START n=%d  np=%d next=%d up=%d mass=%g \n", 
			n,  Tree[n].Npart, Tree[n].DNext, Tree[n].DUp,Tree[n].Mass);
bool tmp = particle_is_inside_node(0, 0);
	add_node(0, 0); // add the first particle by hand
	Tree[0].DUp = -1;
	Tree[0].Bitfield = 0;

	
		

	for (int ipart = 1; ipart < 500; ipart++) {

printf("\n\nIPART=%d NNodes=%zu \n", ipart, NNodes);

		int node = 0;

		for (;;) {


	printf("\nLOOP : n=%d np=%d next=%d level=%d\n", node, Tree[node].Npart, Tree[node].DNext, level(node));

		if (particle_is_inside_node(ipart, node)) { // open
			
	printf("IN : \n");
	
				add_particle_to_node(ipart, node); 
				
				if (Tree[node].Npart == 2) { // refine
				
					int jpart = Tree[node].DNext;
	
	printf("REFINE : node=%d npart=%d jpart=%d\n", 
			node,Tree[node].Npart, jpart);

					Tree[node].DNext = -1;
					
					add_node(jpart, node); // add first daughter
				}  
				
	printf("DECLINE %d->%d  \n", node, node+1);
				node++; // decline

			} else { // skip

	printf("OUT %d : \n", Tree[node].DNext);

				if (Tree[node].DNext < 0 || node == NNodes - 1) {  

	printf("LEAF : \n");

					break; // reached end
				} 

				if (Tree[node].Npart == 1) { // internal leaf

					node++; // might need to split 
	printf("NEXT : +%1\n");
				} else {
			
	printf("NEXT : +%d\n",Tree[node].DNext);
				node += Tree[node].DNext; // skip whole subtree

				}
			}

		} // for (;;)

		if (Tree[node].DNext < 0) // correct internal next node
			Tree[node].DNext = NNodes - node; // only  delta
	
		int parent = node - Tree[node].DUp;

printf("ADD OUT : node=%d npart=%d parent=%d \n", 
			node, Tree[node].Npart, parent);

	add_node(ipart, parent); // add sibling 

	} // for ipart


for (int n=0; n<NNodes; n++) {
printf("n=%d np=%d next=%d up=%d mass=%g level=%d \n", 
n,  Tree[n].Npart, Tree[n].DNext, Tree[n].DUp,Tree[n].Mass, level(n));
print_int_bits64(Tree[n].Bitfield);
}
printf("\n");
for (int ipart = 0; ipart < 500; ipart++) 
print_int_bits64(P[ipart].Peanokey);


exit(0);
	return ;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal 
 */
static inline bool particle_is_inside_node(const int ipart, const int node)
{
	int lvl = level(node);

	if (lvl > 21) // particles closer than 3*21 bit PHKey resolution
		return false; // just add at the end in incoming order (random)

	int shift = 63 - 3 * lvl; // bits 0-5

	uint64_t part_mask = 0x07ULL << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].Peanokey & part_mask) >> shift;

	uint64_t node_triplet = key_fragment(node); 
printf("TEST INSIDE shift=%d level=%d \n", shift, lvl);
print_int_bits64(Tree[node].Bitfield);
print_int_bits64(P[ipart].Peanokey);
print_int_bits64(part_mask);
//print_int_bits64(part_triplet);
//print_int_bits64(node_triplet);

	return (node_triplet == part_triplet);
}
/*
 * We always add nodes at the end of the tree
 */

static inline void add_node(const int ipart, const int parent)
{
	const int node = NNodes++;
	
	Tree[parent].DNext = -1;

	Tree[node].DUp = node - parent;
	Tree[node].DNext = ipart;

	add_particle_to_node(ipart, node); 
	
	uint32_t lvl = level(parent) + 1; // first 5 bits + 1
	
	int shift = 63 - 3*lvl;
printf("ADD ipart=%d n=%d parent=%d lvl=%d shft=%d Par_UP=%d \n", 
		ipart, node, parent,lvl, shift, Tree[parent].DUp);
	uint64_t tmp_key = (P[ipart].Peanokey & (0x07ull << shift));
	uint32_t keyfragment = (uint32_t) (tmp_key >> (shift - 6)); 
	
/*print_int_bits64(P[ipart].Peanokey);
print_int_bits64(0x07ull << shift);
print_int_bits64(keyfragment);
*/
	Tree[node].Bitfield = lvl | keyfragment; 

print_int_bits64(Tree[node].Bitfield);
	return ;
}

static inline void add_particle_to_node(const int ipart, const int node)
{
	Tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	Tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	Tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;
	
	Tree[node].Mass += P[ipart].Mass;

	Tree[node].Npart++;

	/* add yours here */

	return ;
}



/*
 * The level/depth of the node in the tree is saved in 
 * bits 0-5 of the bitmask.
 */

static inline int level(const int node)
{
	return Tree[node].Bitfield & 0x3FULL;
}

/*
 * The 3bit peano key fragment is saved in bits 6-8 of
 * the bitmask.
 */

static inline int key_fragment(const int node)
{
	return (Tree[node].Bitfield & 0x1C0ULL) >> 6;
}

/*
 * The size of the node is just the domain size halfed level times
 */

static inline Float node_size(const int node)
{
	int lvl = level(node);

	return Domain.Size / ((Float) (1UL << lvl));
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
	Tree[0].DUp = -1;
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
