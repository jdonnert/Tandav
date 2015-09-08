#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "domain.h"

static void find_center_of_mass(double *CoM);
static double find_total_mass();

/* 
 * Provide a consistent way of updating/calling different parts 
 * of the code from the main loop. We are already in an OMP parallel
 * environment !
 */

void Update(enum Update_Parameters stage)
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:

		Sig.First_Step = true;

		Set_Current_Cosmology();

		Sim.Total_Mass = find_total_mass();
		find_center_of_mass(Sim.Center_Of_Mass);

		Domain_Decomposition();

		Compute_Acceleration();

		Time_For_Domain_Update();

		Print_Memory_Usage();

		if (Time_For_Snapshot())
			Write_Snapshot();

		break;

	case BEFORE_STEP:

		Write_Logs();

		Sim.Total_Mass = find_total_mass();
		find_center_of_mass(&Sim.Center_Of_Mass);

		break;

	case BEFORE_FIRST_KICK:

		break;

	case BEFORE_SNAPSHOT:

		break;

	case BEFORE_DRIFT:

		break;

	case BEFORE_FORCES:

		Periodic_Constrain_Particles_To_Box(); // PERIODIC

		break;

	case BEFORE_SECOND_KICK:

		break;

	case AFTER_STEP:

		Sig.First_Step = false;

		break;

	default:
		Assert(false, "Update Switch %d not handled", stage);
	}

	return;
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

static double CoM_X = 0, CoM_Y = 0, CoM_Z = 0;

void find_center_of_mass(double CoM_out[3])
{
	#pragma omp single
	CoM_X = CoM_Y = CoM_Z = 0;

	#pragma omp for reduction(+:CoM_X,CoM_Y,CoM_Z,Mass)
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

	for (int i = 0; i < 3; i++)
		CoM_out[i] = global_com[i] / Sim.Center_Of_Mass;

	} // omp single

	return ;

}

