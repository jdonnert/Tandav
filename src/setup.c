#include "globals.h"
#include "timestep.h"

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

