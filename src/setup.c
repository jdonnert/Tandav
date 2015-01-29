#include "globals.h"
#include "timestep.h"

/*
 * Setup extra modules. In particular, allocate static memory blocks here !
 */

void Setup() 
{
	Profile("Setup");

	Setup_Time_Integration();

#ifdef COMOVING
	Setup_Cosmology();

	Setup_Comoving();
#endif // COMOVING
	
	Profile("Setup");

	return ;
}

