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
 * To make the subtree movable, the "next" pointer just gives the difference 
 * between current and next node.
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
	int Up;
} *Tree = NULL;

static inline int level(const int node);
static inline int key_fragment(const int node);
static inline float node_size(const int node);

static size_t NNodes = 0;
static size_t Max_Nodes = 0;
static double Boxsize = 0;

static void add_node(const int ipart, const int node);
static bool part_is_inside_node(const int ipart, const int node);

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
			n,  Tree[n].Npart, Tree[n].DNext, Tree[n].Up,Tree[n].Mass);
bool tmp = part_is_inside_node(0, 0);
	add_node(0, 0); // add the first particle by hand
	Tree[0].Up = -1;
	Tree[0].Bitfield = 0;

	for (int ipart = 1; ipart < 20; ipart++) {

printf("IPART=%d NNodes=%d \n", ipart, NNodes);
		
		Float px = P[ipart].Pos[0] - Domain.Corner[0];
		Float py = P[ipart].Pos[1] - Domain.Corner[1];
		Float pz = P[ipart].Pos[2] - Domain.Corner[2];

		Float mpart = P[ipart].Mass;

		int node = 0;
		int parent = 0;
int cnt=0;
		for (;;) {

if (cnt++ > 10) break;

	printf("LOOP : n=%d np=%d next=%d level=%d\n", node, Tree[node].Npart, Tree[node].DNext, level(node));

		if (part_is_inside_node(ipart, node)) { // open
			
	printf("IN : \n");
				Tree[node].CoM[0] += px * mpart;
				Tree[node].CoM[1] += py * mpart;
				Tree[node].CoM[2] += pz * mpart;
				
				Tree[node].Mass += mpart;
				
				Tree[node].Npart++;

				parent = node;
	
				if (Tree[node].Npart == 2) { // split
				
					int jpart = Tree[node].DNext;
	
	printf("REFINE : node=%d npart=%d jpart=%d\n", node,Tree[node].Npart, jpart);

					Tree[node].DNext = -1;
					
					add_node(jpart, parent); // add daughter
				}  
				
	printf("DECLINE %d->%d  \n", node, node+1);
				node++; // decline

			} else { // skip

				parent = Tree[node].Up;

	printf("OUT %d : \n", Tree[parent].DNext);
				if (Tree[parent].DNext < 0) {
					
	printf("LEAF : \n");

					Tree[node].DNext = NNodes - node; // next give only delta

					break;
				} 
			
				node += Tree[node].DNext;
			}

		} // for (;;)

	
	printf("ADD OUT : node=%d npart=%d \n", node, Tree[node].Npart);
	add_node(ipart, parent); // add sibling or daughter
		

for (int n=0; n<NNodes; n++) {
	printf("n=%d ip=%d np=%d next=%d up=%d mass=%g level=%d \n", 
			n, ipart, Tree[n].Npart, Tree[n].DNext, Tree[n].Up,Tree[n].Mass, level(n));

print_int_bits64(Tree[n].Bitfield);
}

	} // for ipart


exit(0);
	return ;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal 
 */
static bool part_is_inside_node(const int ipart, const int node)
{
	int lvl = level(node);
	int shift = 63 - 3 * lvl; // bits 0-5

	if (shift < 0) // particles closer than 3*21 bit PHKey resolution
		return ipart % 8; // add daughter if true, add sibling if not

	uint64_t part_mask = 0x07ULL << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].Peanokey & part_mask) >> shift;

	uint64_t node_triplet = key_fragment(node); 
printf("TEST INSIDE shift=%d level=%d \n", shift, lvl);
print_int_bits64(Tree[node].Bitfield);
print_int_bits64(part_mask);
print_int_bits64(part_triplet);
print_int_bits64(node_triplet);
print_int_bits64(P[ipart].Peanokey);

	return node_triplet == part_triplet;
}
/*
 * We always add nodes at the end of the tree
 */

static void add_node(const int ipart, const int parent)
{
	const int node = NNodes++;
	
	Tree[parent].DNext = -1;

	Tree[node].Up = parent;
	Tree[node].DNext = ipart;

	Tree[node].CoM[0] = (P[ipart].Pos[0] - Domain.Corner[0]) * P[ipart].Mass;
	Tree[node].CoM[1] = (P[ipart].Pos[1] - Domain.Corner[1]) * P[ipart].Mass;
	Tree[node].CoM[2] = (P[ipart].Pos[2] - Domain.Corner[2]) * P[ipart].Mass;

	Tree[node].Npart = 1;
	Tree[node].Mass = P[ipart].Mass; 

	uint32_t lvl = level(parent) + 1; // first 5 bits + 1
	
	int shift = 63 - 3*lvl;
printf("ADD n=%d parent=%d lvl=%d shft=%d Par_UP=%d \n", 
		node, parent,lvl, shift, Tree[parent].Up);
	uint64_t tmp_key = (P[ipart].Peanokey & (0x07ull << shift));
	uint32_t keyfragment = (uint32_t) (tmp_key >> (shift - 6)); 
	
print_int_bits64(P[ipart].Peanokey);
print_int_bits64(0x07ull << shift);
print_int_bits64(keyfragment);


print_int_bits64(keyfragment);

	Tree[node].Bitfield = lvl | keyfragment; 

print_int_bits64(Tree[node].Bitfield);
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

static inline float node_size(const int node)
{
	int lvl = level(node);

	return Domain.Size[0] / (float) (1UL << lvl);
}


void Init_Tree()
{
	#pragma omp single
	{

	Boxsize = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));
	
	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	if (Tree == NULL)
		Tree = Malloc(Max_Nodes * sizeof(*Tree), "Tree");
	
	Tree[0].CoM[0] = 0; // will hold global CoM of sim
	Tree[0].CoM[1] = 0;
	Tree[0].CoM[2] = 0;

	Tree[0].Npart = 0;
	Tree[0].Up = -1;
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
