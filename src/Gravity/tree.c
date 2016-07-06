#include "tree.h"

#ifdef GRAVITY_TREE

static void check_total_momentum(const bool show_change);

/*
 * Barnes & Hutt gravity tree driver routine
 */

void Gravity_Tree_Acceleration()
{
	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	if (Sig.Prepare_Step) {

		Gravity_Tree_Walk(true); // with BH criterion
	
		Gravity_Tree_Periodic(true); // PERIODIC

		Safe_Last_Accel();
	}

	check_total_momentum(false);

	Gravity_Tree_Walk(false); // relative criterion
		
	Gravity_Tree_Periodic(false); // PERIODIC

	check_total_momentum(true);

	return ;
}
/*
 * Compute total momentum to check the gravity interaction. 
 */

static double Px = 0, Py = 0, Pz = 0, Last = 0;

static void check_total_momentum(const bool show_change)
{
	const int last_DM_part = Task.Npart[0]+Task.Npart[1];

	#pragma omp for reduction(+: Px, Py, Pz)
	for (int ipart = Task.Npart[0]; ipart < last_DM_part; ipart++) {
	
		Px += P.Mass[ipart] * P.Vel[0][ipart];
		Py += P.Mass[ipart] * P.Vel[1][ipart];
		Pz += P.Mass[ipart] * P.Vel[2][ipart];
	}
	
	double ptotal = sqrt( Px*Px + Py*Py + Pz*Pz );

	#pragma omp single
	MPI_Reduce(MPI_IN_PLACE, &ptotal, 1, MPI_DOUBLE, MPI_SUM, MASTER, 
			   MPI_COMM_WORLD);

	double rel_err = (ptotal - Last) / Last;

	if (show_change)
		rprintf("Total change in momentum due to gravity: %g \n", rel_err);

	#pragma omp single nowait
	Last = ptotal;

	return ;
}
#endif
