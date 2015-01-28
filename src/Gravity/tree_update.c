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

	int i = Tree[node].DUp;

	if (D[i].TNode.Level > 0)
		D[i].TNode.Level *= -1; // mark kicked. Will revere after drift

	D[i].TNode.Dp[0] += dp[0] / D[i].TNode.Mass;
	D[i].TNode.Dp[1] += dp[1] / D[i].TNode.Mass;
	D[i].TNode.Dp[2] += dp[2] / D[i].TNode.Mass;

	return ;
}

/*  
 * Advance updated/kicked Treenodes by the system timestep. Then do the same
 * with the top nodes.
 */

static int nUpdate = 0;

void Gravity_Tree_Update_Drift(const double dt)
{
	rprintf("Tree update ");

	#pragma omp for nowait
	for (int i = 0; i < NNodes; i++) {

		if (! Node_Is(UPDATED, i))
			continue;

		Tree[i].CoM[0] += dt * Tree[i].Dp[0];
		Tree[i].CoM[1] += dt * Tree[i].Dp[1];
		Tree[i].CoM[2] += dt * Tree[i].Dp[2];

		Tree[i].Dp[0] = Tree[i].Dp[1] = Tree[i].Dp[2] = 0;

		Node_Clear(UPDATED, i);
	} 

	nUpdate = 0;

	#pragma omp for reduction(+:nUpdate)
	for (int i = 0; i < NTop_Nodes; i++) {
	
		if (D[i].TNode.Level > 0)
			continue;

		D[i].TNode.CoM[0] += dt / D[i].TNode.Dp[0];
		D[i].TNode.CoM[1] += dt / D[i].TNode.Dp[1];
		D[i].TNode.CoM[2] += dt / D[i].TNode.Dp[2];

		D[i].TNode.Level *= -1; // reverse "updated" flag

		nUpdate++;
	}

	rprintf("Done. Moved %d top nodes \n", nUpdate);

	return ;
}


