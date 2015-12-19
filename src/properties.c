#include "globals.h"

static void find_center_of_mass(double CoM_out[3]);
static void find_total_mass(double mass_out[1]);

void Compute_Global_Simulation_Properties()
{
	find_total_mass(&Sim.Total_Mass);

	find_center_of_mass(&Sim.Center_Of_Mass[0]);

	return ;
}

static double Total_Mass = 0;

static void find_total_mass(double mass_out[1])
{
	#pragma omp single
	Total_Mass = 0;

	#pragma omp for reduction(+:Total_Mass)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		Total_Mass += P.Mass[ipart];

	#pragma omp single
	{

	MPI_Allreduce(MPI_IN_PLACE, &Total_Mass, 1, MPI_DOUBLE, MPI_SUM,
				  MPI_COMM_WORLD);

	mass_out[0] = Total_Mass;

	} // omp single

	return ;
}

static double CoM_X = 0, CoM_Y = 0, CoM_Z = 0; // can't reduce on array

static void find_center_of_mass(double CoM_out[3])
{
	#pragma omp single
	CoM_X = CoM_Y = CoM_Z = 0;

	#pragma omp for reduction(+:CoM_X,CoM_Y,CoM_Z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		CoM_X += P.Mass[ipart] * P.Pos[0][ipart];
		CoM_Y += P.Mass[ipart] * P.Pos[1][ipart];
		CoM_Z += P.Mass[ipart] * P.Pos[2][ipart];
	}

	#pragma omp single
	{

	double global_com[3] = { CoM_X, CoM_Y, CoM_Z  };

	MPI_Allreduce(MPI_IN_PLACE, &global_com, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	CoM_out[0] = global_com[0] / Sim.Total_Mass;
	CoM_out[1] = global_com[1] / Sim.Total_Mass;
	CoM_out[2] = global_com[2] / Sim.Total_Mass;

	} // omp single

	return ;

}

