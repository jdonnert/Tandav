#include "globals.h"
#include "proto.h"
#include "update.h"

void Update(enum Update_Parameters stage) 
{

	if (stage == BEFORE_MAIN_LOOP) {
		
		Print_Memory_Usage();
	
	}

	if (stage == AFTER_FIRST_KICK) {

	}

	if (stage == AFTER_DRIFT) {

	}

	if (stage == AFTER_SECOND_KICK) {

	}

	return;
}
