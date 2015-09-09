#include "globals.h"

static void find_center_of_mass(double CoM_out[3]);
static double find_total_mass();

void Compute_Global_Simulation_Properties()
{
	Sim.Total_Mass = find_total_mass();

	find_center_of_mass(Sim.Center_Of_Mass);

	return ;
}

static double total_mass = 0;

static double find_total_mass()
{
	#pragma omp single
	total_mass = 0;

	#pragma omp for reduction(+:total_mass)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		total_mass += P[ipart].Mass;

	#pragma omp single
	MPI_Allreduce(MPI_IN_PLACE, &total_mass, 1, MPI_DOUBLE, MPI_SUM,
				  MPI_COMM_WORLD);

	return total_mass;
}

static double CoM_X = 0, CoM_Y = 0, CoM_Z = 0; // can't reduce on array

static void find_center_of_mass(double CoM_out[3])
{
	#pragma omp single
	CoM_X = CoM_Y = CoM_Z = 0;

	#pragma omp for reduction(+:CoM_X,CoM_Y,CoM_Z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		CoM_X += P[ipart].Mass * P[ipart].Pos[0];
		CoM_Y += P[ipart].Mass * P[ipart].Pos[1];
		CoM_Z += P[ipart].Mass * P[ipart].Pos[2];
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

