#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "io/io.h"
#include "peano.h"

static void find_global_boxsize();
		
/* provide a consistent way of updating/calling different parts 
 * of the code from the main loop without cluttering */

void Update(enum Update_Parameters stage) 
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:

#ifndef PERIODIC
		find_global_boxsize();
#endif
		Sort_Particles_By_Peano_Key();
		
		Compute_Acceleration();
		
		Print_Memory_Usage();

		#pragma omp barrier

		if (Time.Begin == Time.Next_Snap) { 

			Write_Snapshot();

			#pragma omp single
			Time.Next_Snap += Time.Bet_Snap;
		}
		
		break;

	case BEFORE_STEP:
		
		Write_Logs();

#ifndef PERIODIC
		find_global_boxsize();
#endif

		break;

	case BEFORE_SNAPSHOT:

		break;
		
	case BEFORE_DRIFT:
		
		break;

	case BEFORE_FORCES:

		Sort_Particles_By_Peano_Key();

		break;

	case BEFORE_SECOND_KICK:

		break;
		
	default:
		Assert(false, "Update Stage %d not handled", stage);
	}

	return;
}


static void find_global_boxsize()
{
	double local_max[3] = { 0 };

	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		local_max[0] = fmax(local_max[0], fabs(P[ipart].Pos[0]));
		local_max[1] = fmax(local_max[1], fabs(P[ipart].Pos[1]));
		local_max[2] = fmax(local_max[2], fabs(P[ipart].Pos[2]));
	}

	#pragma omp single 
	{

	MPI_Allreduce(&local_max, &Sim.Boxsize, 3, MPI_DOUBLE, MPI_MAX,
			MPI_COMM_WORLD);

	Sim.Boxsize[0] *= 2;
	Sim.Boxsize[1] *= 2;
	Sim.Boxsize[2] *= 2;
	
	}

	return ;
}
