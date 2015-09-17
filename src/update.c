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

		Set_Current_Cosmology();

		Compute_Global_Simulation_Properties();

		Domain_Decomposition();

		Compute_Acceleration();

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


