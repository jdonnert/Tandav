#include "globals.h"
#include "domain.h"
#include "timestep.h"
#include "properties.h"

/*
 * Setup extra modules. In particular, allocate static memory blocks here !
 */

void Setup()
{
	Profile("Setup");

	Setup_Cosmology(); // COMOVING

	Setup_Time_Integration();

	Setup_Comoving(); // COMOVING
		
	Compute_Global_Simulation_Properties();

	Setup_Domain_Decomposition();

	Print_Memory_Usage();

	Profile("Setup");

	return ;
}

