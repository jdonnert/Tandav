#include "globals.h"
#include "tree.h"
#include "domain.h"

#define NODES_PER_PARTICLE 0.7 // Springel 2005

struct Tree_Node {
	Float Center[3];
	Float Size;
	int Npart;
	float Mass;
	int Next;		
	int Up;
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

/*	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		Float px = P[ipart].Pos[0] - Domain_Corner[0];
		Float py = P[ipart].Pos[1] - Domain_Corner[1];
		Float pz = P[ipart].Pos[2] - Domain_Corner[2];

		int node = 0;

		for (;;) {
			
			if (is_inside(px, py, pz, node)) {
			
				Tree[node].Mass += P[ipart].Mass;
	
				if (Tree[node].Npart++ < 2) // refine
					break; 
				
				node++; // decline

			} else { // skip

				if (Tree[node].Next < 0) {
					
					node = NNodes; // add at the end
					
					break;
				} 
			
				node += Tree[node].Next;
			}

		} // for (;;)

		if (Tree[node].Npart == 2) { // make space
			
			int jpart = Tree[node].Next;
					
			add_node(jpart, node++, node);
		}
		add_node(ipart, node++, node);
		
	} // for ipart
*/
	return ;
}

static bool is_inside(const Float px, const Float py, const Float pz,
		const int node)
{
	const Float ds = Tree[node].Size * 0.5;

	const Float dx = fabs(px - Tree[node].Center[0]);
	const Float dy = fabs(py - Tree[node].Center[1]);
	const Float dz = fabs(pz - Tree[node].Center[2]);

	return (dx <= ds) && (dy <= ds) && (dz <= ds) ;
}

/* 
 * The node center is determined by the location of the occupying particle.
 */

static void add_node(const int ipart, const int parent, const int node)
{

	Tree[node].Up = parent;
	Tree[node].Next = ipart;

	Tree[node].Size = Tree[parent].Size * 0.5;

	Float px = P[ipart].Pos[0] - Domain_Corner[0];
	Float py = P[ipart].Pos[1] - Domain_Corner[1];
	Float pz = P[ipart].Pos[2] - Domain_Corner[2];

	Float dx = Sign(px - Tree[parent].Center[0]) * 0.5 * Tree[node].Size;
	Float dy = Sign(py - Tree[parent].Center[1]) * 0.5 * Tree[node].Size;
	Float dz = Sign(pz - Tree[parent].Center[2]) * 0.5 * Tree[node].Size;

	Tree[node].Center[0] = Tree[parent].Center[0] + dx;
	Tree[node].Center[1] = Tree[parent].Center[1] + dy;
	Tree[node].Center[2] = Tree[parent].Center[2] + dz;

	Tree[node].Npart = 1;
	Tree[node].Mass = P[ipart].Mass; 

	NNodes++;

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


