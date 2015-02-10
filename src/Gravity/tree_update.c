#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

/*
 * Dynamically update the tree with the kicks of this timestep. We
 * start at the parent of the current particle and walk the tree backwards
 * until we reach the top node. Then we kick the top node itself, and mark it
 * as kicked by multiplying the level with -1
 */

void Gravity_Tree_Update_Kicks(const int ipart, const double dt)
{
	Float m_dt = P[ipart].Mass * dt; // kick tree nodes

	const Float dp[3] = { m_dt*P[ipart].Acc[0],
						  m_dt*P[ipart].Acc[1], 
						  m_dt*P[ipart].Acc[2] };
	
	int i = 0;
	int node = P[ipart].Tree_Parent;

//printf("FK Update ipart=%d parent %d \n", ipart, node);

	if (node >= 0) {
	
		while (!Node_Is(TOP, node)) { // topnode & tree
//printf("node %d next %d up %d  npart %d level %d \n", node, Tree[node].DNext, Tree[node].DUp, Tree[node].Npart, Level(node));
	
			#pragma omp atomic
			Tree[node].Dp[0] += dp[0] / Tree[node].Mass;
			#pragma omp atomic
			Tree[node].Dp[1] += dp[1] / Tree[node].Mass;
			#pragma omp atomic
			Tree[node].Dp[2] += dp[2] / Tree[node].Mass;
		
			#pragma omp atomic
			Tree[node].Bitfield |= 1UL << UPDATED; // = Node_Set(UPDATED,node)

			node -= Tree[node].DUp;
		} // while not top node
//printf("node %d next %d up %d  npart %d level %d \n", node, Tree[node].DNext, Tree[node].DUp, Tree[node].Npart, Level(node));

		i = Tree[node].DUp;

	} else { // topnode only

		i = -node + 1;
		//printf("Top node only ipart %d node %d \n", ipart, i);
	}
	
	if (D[i].TNode.Level > 0)
		D[i].TNode.Level *= -1; // mark kicked. Will reverse after drift

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
	rprintf("Tree update. ");

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

		D[i].TNode.CoM[0] += dt * D[i].TNode.Dp[0];
		D[i].TNode.CoM[1] += dt * D[i].TNode.Dp[1];
		D[i].TNode.CoM[2] += dt * D[i].TNode.Dp[2];

		D[i].TNode.Dp[0] = D[i].TNode.Dp[1] = D[i].TNode.Dp[2] = 0;
			
		D[i].TNode.Level *= -1; // reverse "updated" flag

		nUpdate++;
	}

	rprintf("Moved %d top nodes \n", nUpdate);

	return ;
}


