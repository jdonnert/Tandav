#include "globals.h"
#include "domain.h"
#include "timestep.h"
#include "properties.h"
#include "Gravity/gravity.h"

static void sanity_check_simulation_setup();

/*
 * Setup extra modules. In particular, allocate static memory blocks here !
 */

void Setup()
{
	Profile("Setup");
		
	Setup_Time_Integration();

	Setup_Comoving(); // COMOVING

	Compute_Global_Simulation_Properties();

	Setup_Domain_Decomposition();
	
	Setup_Gravity_Tree(); // GRAVITY_TREE

	/* Add yours above */

	sanity_check_simulation_setup();
	
	Print_Memory_Usage();

	Profile("Setup");

	return ;
}


static void sanity_check_simulation_setup()
{

#ifdef COMOVING

	double omega_matter = Sim.Total_Mass / p3(Sim.Boxsize[0]) 
							/ Cosmo.Critical_Density;

	Warn(fabs(omega_matter-Cosmo.Omega_Matter) > 1e-3, 
			"Matter content of box is inconsistent with Omega_Matter in "
			"parameter file. \n  Run %g, ParFile %g ", 
			omega_matter, Cosmo.Omega_Matter );

#endif // COMOVING

	return ;
}
