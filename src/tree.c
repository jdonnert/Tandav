#include "globals.h"
#include "tree.h"
#include "domain.h"

#define NODES_PER_PARTICLE 0.7 // Springel 2005

struct Tree_Node {
	Float Center[3];
	Float Size;
	int Npart;
	int Next;		
	int Up;
	float Mass;
} *Tree = NULL;

static size_t NNodes = 0;
static size_t Max_Nodes = 0;
static double Boxsize = 0;

void Build_Tree()
{
	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		Float px = P[ipart].Pos[0] - Domain_Corner[0];
		Float py = P[ipart].Pos[1] - Domain_Corner[1];
		Float pz = P[ipart].Pos[2] - Domain_Corner[2];

		int node = 0;
		int parent = -1;
		bool part_done = false;

		for (;;) {
			
			if (is_inside(px, py, pz, node)) {
			
				Tree[node].Npart++;
				Tree[node].Mass += P[ipart].Mass;
	
				if (Tree[node].Npart > 2) // decline
					node += 1;
				else  // refine
					break;

				continue;

			} else { // skip

				node += Tree[node].Next;
			}
		
		} // while

		Tree[node].size = Tree[parent].Size * 0.5;

		Tree[node].Up = parent;
		Tree[node].Center[0] = Tree[parent].Center[0]+sign*Tree[node].size;
		Tree[node].Center[1] = Tree[parent].Center[1]+sign*Tree[node].size;
		Tree[node].Center[2] = Tree[parent].Center[2]+sign*Tree[node].size;

		Tree[node].Npart = 1;
		Tree[node].Mass = P[ipart].Mass; 
	}

	return ;
}

void Init_Tree()
{
#pragma omp single
	{

	Boxsize = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));
	
	Max_Nodes = Task.Npart_Total_Max*NODES_PER_PARTICLE;
	
	Tree = Malloc(Max_Nodes * sizeof(*Tree), "Tree")
	
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


