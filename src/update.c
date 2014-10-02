#include "globals.h"
#include "update.h"
#include "accel.h"
#include "timestep.h"
#include "io/io.h"

/* provide a consistent way of updating/calling different parts 
 * of the code from the main loop without cluttering */

void Update(enum Update_Parameters stage) 
{
	switch (stage) {

	case BEFORE_MAIN_LOOP:
		
		Compute_Acceleration();
		
		Print_Memory_Usage();

		if (Time.Begin == Time.Next_Snap) { 

			Write_Snapshot();

			Time.Next_Snap += Time.Bet_Snap;
		}
		
		break;

	case BEFORE_FIRST_KICK:
	
		break;

	case AFTER_FIRST_KICK:
		
		break;

	case BEFORE_SNAPSHOT:

		break;
		
	case BEFORE_FORCES:
		
		break;

	case BEFORE_SECOND_KICK:

		break;
		
	case AFTER_SECOND_KICK:

		Write_Logs();

		break;
		
	default:
		Assert(0, "Update Stage %d not handled", stage);
	}

	return;
}


