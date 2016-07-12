#include "update.h"

/* 
 * Provide a consistent way of updating/calling different parts 
 * of the code from the main loop. We are already in an OMP parallel
 * environment !
 */

void Update(enum Update_Parameters stage)
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:

		Sig.Prepare_Step = true;

		Set_Current_Cosmology(Time.Current); // COMOVING

		Domain_Decomposition();

		Compute_Acceleration();

		Time_For_Domain_Update();

		Print_Memory_Usage();

		if (Time_For_Snapshot())
			IO_Write_Snapshot();

		Sig.Prepare_Step = false;

		Sig.First_Step = true;
		
		break;

	case BEFORE_STEP:
		
		Compute_Current_Simulation_Properties();

		Write_Logs();

		break;

	case BEFORE_FIRST_KICK:

		break;

	case BEFORE_SNAPSHOT:

		break;

	case BEFORE_DRIFT:

		break;

	case AFTER_DRIFT:

		Compute_Current_Simulation_Properties();

		Periodic_Constrain_Particles_To_Box(); // PERIODIC

		break;

	case BEFORE_DOMAIN_UPDATE: 

		Gravity_Tree_Free(); // GRAVITY_TREE
		
		Gravity_FMM_Free(FMM); // GRAVITY_FMM
	
		break;

	case BEFORE_FORCES:

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


// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
