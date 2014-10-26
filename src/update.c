#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "IO/io.h"
#include "domain.h"
#include "Gravity/gravity.h"

/* 
 * Provide a consistent way of updating/calling different parts 
 * of the code from the main loop
 */

void Update(enum Update_Parameters stage) 
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:

#pragma omp parallel
		{

		Update(FORCES);
		
		Print_Memory_Usage();

		Compute_Acceleration();

		} // omp parallel
		
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

	case BEFORE_SNAPSHOT:

		break;
		
	case BEFORE_DRIFT:
		
		break;

	case FORCES:

		Domain_Decomposition();

#ifdef GRAVITY_TREE
		Build_Tree();
#endif

		break;

	case BEFORE_SECOND_KICK:

		break;
		
	default:
		Assert(false, "Update Stage %d not handled", stage);
	}

	return;
}

