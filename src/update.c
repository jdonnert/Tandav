#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "domain.h"

/* 
 * Provide a consistent way of updating/calling different parts 
 * of the code from the main loop. We are already in an OMP parallel
 * environment !
 */

void Update(enum Update_Parameters stage) 
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:
		
#ifdef COMOVING
		Set_Current_Cosmology();
#endif

		Domain_Decomposition();


		Print_Memory_Usage();

		Compute_Acceleration();

		Print_Memory_Usage();

		if (Time.Begin == Time.Next_Snap) { 
			
			Write_Snapshot();
			
			#pragma omp single
			Time.Next_Snap += Time.Bet_Snap;
		}
		
		Print_Memory_Usage();

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

	case BEFORE_FORCES:

#ifdef COMOVING
		Set_Current_Cosmology();
#endif
		break;

	case BEFORE_SECOND_KICK:

		break;
	
	case AFTER_STEP:

		#pragma omp single
		Sig.First_Step = false;

		break;

	default:
		Assert(false, "Update Switch %d not handled", stage);
	}

	return;
}

