/*
 * This builds the tree 
 */

#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#define NODES_PER_PARTICLE 0.7 // Springel 2005

struct Tree_Node {
	Float CoM[3];
	uint32_t Bitmask; // 5b:level, 3b:key, 24b:free
	int Npart;
	float Mass;
	int Next;		
	int Up;
} *Tree = NULL;

static size_t NNodes = 0;
static size_t Max_Nodes = 0;
static double Boxsize = 0;

static void add_node(const int, const int);
static bool is_inside(const Float, const Float, const Float, const int);

void Build_Tree()
{
	memset(Tree, 0, Max_Nodes * sizeof(*Tree));

	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		Float px = P[ipart].Pos[0] - Domain.Corner[0];
		Float py = P[ipart].Pos[1] - Domain.Corner[1];
		Float pz = P[ipart].Pos[2] - Domain.Corner[2];

		Float mpart = P[ipart].Mass;

		int node = 0;
		int parent = 0;

		for (;;) {
			
			if (is_inside(px, py, pz, node)) { // open
			
				Tree[node].CoM[0] += px * mpart;
				Tree[node].CoM[1] += py * mpart;
				Tree[node].CoM[2] += pz * mpart;
				
				Tree[node].Mass += mpart;
				
				Tree[node].Npart++;

				parent = node;
	
				if (Tree[node].Npart == 2) // split
					break; 
				
				node++; // decline

			} else { // skip

				if (Tree[node].Next < 0) {
					
					parent = Tree[node].Up;

					break;
				} 
			
				node += Tree[node].Next;
			}

		} // for (;;)

		if (Tree[node].Npart == 2) { // make space
			
			int jpart = Tree[node].Next;
					
			add_node(jpart, parent);
		}

		add_node(ipart, parent);
		
	} // for ipart

	return ;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal 
 */
static bool is_inside(const int ipart, const int node)
{
	int shift = 3*(Tree[node].Bitmask & 0x1F); 
	uint64_t part_mask = 0x07 << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].PeanoKey & part_mask) >> shift;

	uint64_t node_triplet = Tree[node].Bitmask & 0xE0; 

	return node_triplet == part_triplet;
}

/*
 * We always add nodes at the end
 */
static void add_node(const int ipart, const int parent)
{
	const int node = NNodes++;

	Tree[node].Up = parent;
	Tree[node].Next = ipart;

	Tree[node].CoM[0] = (P[ipart].Pos[0] - Domain.Corner[0]) * P[ipart].Mass;
	Tree[node].CoM[1] = (P[ipart].Pos[1] - Domain.Corner[1]) * P[ipart].Mass;
	Tree[node].CoM[2] = (P[ipart].Pos[2] - Domain.Corner[2]) * P[ipart].Mass;

	Tree[node].Npart = 1;
	Tree[node].Mass = P[ipart].Mass; 

	uint32_t level = (Tree[parent].Bitmask & 0x1F) + 1; // first 5 bits + 1
	
	uint32_t keyfragment = P[ipart].PeanoKey & (0x1F << 3*level);
	
	keyfragment >>= 3*level - 5; // leave space for level

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
	
	Tree[0].Center[0] = 0.5 * Boxsize;
	Tree[0].Center[1] = 0.5 * Boxsize;
	Tree[0].Center[2] = 0.5 * Boxsize;

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
void Accel_Gravity_Tree(const int ipart, double acc*, double *pot)
{

	return ;
}
