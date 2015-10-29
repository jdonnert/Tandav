#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "properties.h"

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

		Set_Current_Cosmology(Time.Current); // COMOVING

		Compute_Global_Simulation_Properties();

		Domain_Decomposition();

		Compute_Acceleration();

#pragma omp single
for (int ipart=0; ipart<Sim.Npart_Total; ipart++)
	if (P[ipart].ID == 1)
		printf("ACC UPDATE DONE ipart=%d, ID=%d pos=%g %g %g, vel=%g %g %g acc= %g %g %g \n",ipart, P[ipart].ID, P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2], P[ipart].Vel[0],P[ipart].Vel[1],P[ipart].Vel[2],P[ipart].Grav_Acc[0] ,P[ipart].Grav_Acc[1],P[ipart].Grav_Acc[2]);

//exit(0);
		Time_For_Domain_Update();

		Print_Memory_Usage();

		if (Time_For_Snapshot())
			Write_Snapshot();

		break;

	case BEFORE_STEP:

		Write_Logs();

		break;

	case BEFORE_FIRST_KICK:

		break;

	case BEFORE_SNAPSHOT:

		break;

	case BEFORE_DRIFT:

		break;

	case BEFORE_DOMAIN:

		Compute_Global_Simulation_Properties();

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
		Assert(false, "Update stage %d not handled", stage);
	}

	return;
}


