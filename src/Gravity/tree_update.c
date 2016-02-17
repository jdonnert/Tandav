#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "tree.h"

#ifdef GRAVITY_TREE

/*
 * Dynamically update the tree with the kicks of this timestep. We
 * start at the parent of the current particle and walk the tree backwards
 * until we reach the top node. Then we kick the top node itself, and mark it
 * as kicked by multiplying the level with -1
 */

void Gravity_Tree_Update_Kicks(const int ipart, const double dt)
{
	Float m_dt = P.Mass[ipart] * dt; // kick tree nodes

	const Float dp[3] = { m_dt*P.Acc[0][ipart],
						  m_dt*P.Acc[1][ipart],
						  m_dt*P.Acc[2][ipart] };
	int i = 0;
	int node = P.Tree_Parent[ipart];

	if (node >= 0) { // kick sub tree nodes

		while (! Node_Is(TOP, node)) {

			#pragma omp atomic update
			Tree[node].Dp[0] += dp[0] / Tree[node].Mass;
			#pragma omp atomic update
			Tree[node].Dp[1] += dp[1] / Tree[node].Mass;
			#pragma omp atomic update
			Tree[node].Dp[2] += dp[2] / Tree[node].Mass;

			#pragma omp atomic update
			Tree[node].Bitfield |= 1UL << UPDATED; // = Node_Set(UPDATED,node)

			node -= Tree[node].DUp;
		} // while 

		i = Tree[node].DUp;

	} else  // kick top node only
		i = -node - 1;

	#pragma omp critical
	if (D[i].TNode.Level > 0)
		D[i].TNode.Level *= -1; // mark kicked. Will reverse after drift

	#pragma omp atomic update
	D[i].TNode.Dp[0] += dp[0] / D[i].TNode.Mass;
	#pragma omp atomic update
	D[i].TNode.Dp[1] += dp[1] / D[i].TNode.Mass;
	#pragma omp atomic update
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
	#pragma omp single
	nUpdate = 0;

	#pragma omp for nowait
	for (int i = 0; i < NNodes; i++) {

		if (Node_Is(UPDATED, i)) {

			Tree[i].CoM[0] += dt * Tree[i].Dp[0];
			Tree[i].CoM[1] += dt * Tree[i].Dp[1];
			Tree[i].CoM[2] += dt * Tree[i].Dp[2];

			Tree[i].Dp[0] = Tree[i].Dp[1] = Tree[i].Dp[2] = 0;

			Node_Clear(UPDATED, i);
		}
	}

	#pragma omp for simd reduction(+:nUpdate)
	for (int i = 0; i < NTop_Nodes; i++) {

		if (D[i].TNode.Level < 0) {

			D[i].TNode.CoM[0] += dt * D[i].TNode.Dp[0];
			D[i].TNode.CoM[1] += dt * D[i].TNode.Dp[1];
			D[i].TNode.CoM[2] += dt * D[i].TNode.Dp[2];

			D[i].TNode.Dp[0] = D[i].TNode.Dp[1] = D[i].TNode.Dp[2] = 0;

			D[i].TNode.Level *= -1; // reverse "updated" flag

			nUpdate++;
		}
	}

	rprintf("Tree update: Moved %d top nodes \n", nUpdate);

	return ;
}

#endif // GRAVITY_TREE


