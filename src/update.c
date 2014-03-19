#include "globals.h"
#include "proto.h"
#include "update.h"

/* provide a consistent way of updating/calling different parts of the code
 * This file should contain only function calls */
void Update(enum Update_Parameters stage) 
{
	switch (stage) {

		case BEFORE_MAIN_LOOP:
		
			Print_Memory_Usage();
		
			break;

		case AFTER_FIRST_KICK:
		
			break;
		
		case AFTER_DRIFT:
		
			break;
		
		case AFTER_SECOND_KICK:

			break;
		
		default:
			Assert(0, "Update Stage %d not handled", stage);
	}

	return;
}
