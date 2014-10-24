/*
 * This builds the tree 
 */

#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#define NODES_PER_PARTICLE 0.7 // Springel 2005

struct Tree_Node {
	Float CoM[3];
	uint32_t Bitmask; // bit 0-4:level, 5-7:key, 8-32:free
	int Npart;
	float Mass;
	int Next;		
	int Up;
} *Tree = NULL;

static size_t NNodes = 0;
static size_t Max_Nodes = 0;
static double Boxsize = 0;

static void add_node(const int ipart, const int node);
static bool is_inside(const int ipart, const int node);

static void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--)
		printf("%llu", (val & (1ULL << i) ) >> i);
	
	printf("\n");fflush(stdout);

	return ;
}


void Build_Tree()
{
	memset(Tree, 0, Max_Nodes * sizeof(*Tree));

	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/

for (int n=0; n<NNodes; n++)
	printf("START n=%d  np=%d next=%d up=%d mass=%g \n", 
			n,  Tree[n].Npart, Tree[n].Next, Tree[n].Up,Tree[n].Mass);

	for (int ipart = 0; ipart < 2; ipart++) {
	printf("IPART=%d NNodes=%d \n", ipart, NNodes);
		Float px = P[ipart].Pos[0] - Domain.Corner[0];
		Float py = P[ipart].Pos[1] - Domain.Corner[1];
		Float pz = P[ipart].Pos[2] - Domain.Corner[2];

		Float mpart = P[ipart].Mass;

		int node = 0;
		int parent = 0;

		for (;;) {
	printf("LOOP : n=%d \n", node);
			if (is_inside(ipart, node)) { // open
			
	printf("IN : \n");
				Tree[node].CoM[0] += px * mpart;
				Tree[node].CoM[1] += py * mpart;
				Tree[node].CoM[2] += pz * mpart;
				
				Tree[node].Mass += mpart;
				
				Tree[node].Npart++;

				parent = node;
	
				if (Tree[node].Npart <= 2) // split
					break; 
				
	printf("DECLINE : \n");
				node++; // decline

			} else { // skip

	printf("OUT : \n");
				if (Tree[node].Next < 0) {
					
	printf("LEAF : \n");
					parent = Tree[node].Up;

					Tree[node].Next = NNodes - node;

					break;
				} 
			
				node += Tree[node].Next;
			}

		} // for (;;)

		if (Tree[node].Npart == 2) { // make space
			
	printf("REFINE : node=%d npart=%d \n", node, Tree[node].Npart);
			int jpart = Tree[node].Next;
					
			add_node(jpart, parent);

		} 
	
	printf("ADD : node=%d npart=%d \n", node, Tree[node].Npart);
		add_node(ipart, parent);
		

for (int n=0; n<NNodes; n++)
	printf("n=%d ip=%d np=%d next=%d up=%d mass=%g \n", 
			n, ipart, Tree[n].Npart, Tree[n].Next, Tree[n].Up,Tree[n].Mass);
	
	} // for ipart


exit(0);
	return ;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal 
 */
static bool is_inside(const int ipart, const int node)
{
	int shift = 63 - 3 * ((Tree[node].Bitmask & 0x1F) >> 6); // extract bits 6-8

	uint64_t part_mask = 0x07 << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].Peanokey & part_mask) >> shift;

	uint64_t node_triplet = Tree[node].Bitmask & 0xE0; 

print_int_bits64(Tree[node].Bitmask);
print_int_bits64(part_mask);
print_int_bits64(part_triplet);
print_int_bits64(node_triplet);
print_int_bits64(P[ipart].Peanokey);

	return node_triplet == part_triplet;
}

/*
 * We always add nodes at the end
 */
static void add_node(const int ipart, const int parent)
{
	const int node = NNodes++;

printf("ADD n=%d parent=%d \n", node, parent);

	Tree[node].Up = parent;
	Tree[node].Next = ipart;

	Tree[node].CoM[0] = (P[ipart].Pos[0] - Domain.Corner[0]) * P[ipart].Mass;
	Tree[node].CoM[1] = (P[ipart].Pos[1] - Domain.Corner[1]) * P[ipart].Mass;
	Tree[node].CoM[2] = (P[ipart].Pos[2] - Domain.Corner[2]) * P[ipart].Mass;

	Tree[node].Npart = 1;
	Tree[node].Mass = P[ipart].Mass; 

	uint32_t level = (Tree[parent].Bitmask & 0x1F) + 1; // first 5 bits + 1
	
	uint32_t keyfragment = P[ipart].Peanokey & (0x1F << 63 - 3*level);
	
	keyfragment >>= 3*level - 5; // leave space for level

print_int_bits64(Tree[node].Bitmask);

	Tree[node].Bitmask = level | keyfragment; 

	return ;
}

void Init_Tree()
{
	#pragma omp single
	{

	Boxsize = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));
	
	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	Tree = Malloc(Max_Nodes * sizeof(*Tree), "Tree");
	
	Tree[0].CoM[0] = 0; // will hold global CoM of sym
	Tree[0].CoM[1] = 0;
	Tree[0].CoM[2] = 0;

	Tree[0].Npart = 0;
	Tree[0].Up = -1;
	Tree[0].Mass = 0;
	Tree[0].Next = -1;

	NNodes = 1;

	} // omp single

	return ;
}

/*
 * This function computes the gravitational acceleration by walking the tree  
 */
void Accel_Gravity_Tree(const int ipart, double *acc, double *pot)
{

	return ;
}
