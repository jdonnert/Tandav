#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

/*
 * Dynamically update the tree with the kicks of this timestep. We
 * start at the parent of the current particle and walk the tree backwards
 * until we reach the first top node, which corresponds to the bunchleave of
 * this subtree.
 */

void Gravity_Tree_Update_Kicks(const Float dp[3], const int parent)
{
	int node = parent;
	
	do {

		#pragma omp critical
		{
		Tree[node].Dp[0] += dp[0] / Tree[node].Mass;
		Tree[node].Dp[1] += dp[1] / Tree[node].Mass;
		Tree[node].Dp[2] += dp[2] / Tree[node].Mass;
		
		Node_Set(UPDATED, node);
		}

		node -= Tree[node].DUp;

	} while (!Node_Is(TOP, node));

	return ;
}

/*  
 * Advance the Treenodes by the system timestep, even though we just drift
 * the 
 */

void Gravity_Tree_Update_Drift(const double dt)
{

	rprintf("Tree update \n");

	#pragma omp for
	for (int i = 0; i < NNodes; i++) {

		if (! Node_Is(UPDATED, i))
			continue;

		Tree[i].CoM[0] += dt * Tree[i].Dp[0];
		Tree[i].CoM[1] += dt * Tree[i].Dp[1];
		Tree[i].CoM[2] += dt * Tree[i].Dp[2];

		Tree[i].Dp[0] = Tree[i].Dp[1] = Tree[i].Dp[2] = 0;

		Node_Clear(UPDATED, i);
	} 

	return ;
}

/*
 * Communicate the kicks of the bunch leaves and kick the topnode tree.
 * We find all changed bunch leave momenta
 */

void Gravity_Tree_Update_Topnode_Kicks()
{
	// Find updated bunch leaves
	// MPI allscatter bunch leave momenta
	// update topnode tree
	
	return ;
}

