#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

/*
 * Dynamically update the tree with the kicks of this timestep
 */

void Gravity_Tree_Update_Kicks(const Float dv[3], const int parent)
{
	return ;
	int node = parent;

	while (!Node_Is(TOP, node)) {
		
		
	}

	// add to bunchleave kicks


	return ;
}

/*  
 * Advance the Treenodes by the system timestep
 */

void Gravity_Tree_Update_Drift(const double dt)
{
return ;
	#pragma omp for
	for (int i = 0; i < NNodes; i++) {

		Tree[i].CoM[0] += dt * Tree[i].Vel_CoM[0];
		Tree[i].CoM[1] += dt * Tree[i].Vel_CoM[1];
		Tree[i].CoM[2] += dt * Tree[i].Vel_CoM[2];
	}

	return ;
}

/*
 * Communicate the kicks of the bunch leaves and kick the topnode tree
 */

void Gravity_Tree_Update_Topnode_Kicks()
{
	// zero bunch leave kicks
	// MPI reduce bunch leave kicks
	


	return ;
}

