#include "globals.h"
#include "domain.h"
#include "timestep.h"
#include "properties.h"
#include "Gravity/gravity.h"

static void sanity_check_simulation();

/*
 * Setup extra modules. In particular, allocate static memory blocks here !
 */

void Setup()
{
	Profile("Setup");

	Setup_Periodic(); // PERIODIC

	Setup_Cosmology(); // COMOVING

	Setup_Time_Integration();

	Setup_Comoving(); // COMOVING
		
	Print_Memory_Usage();
	
	Compute_Global_Simulation_Properties();

	Setup_Domain_Decomposition();
	
	Setup_Gravity_Tree(); // GRAVITY_TREE

	Print_Memory_Usage();

	sanity_check_simulation();

	Profile("Setup");

	return ;
}


static void sanity_check_simulation()
{

#ifdef COMOVING

	double omega_matter = Sim.Total_Mass / p3(Sim.Boxsize[0]) 
		/ (3 * p2(Cosmo.Hubble_Parameter) / (8 * Pi * Const.Gravity) );

	Warn( fabs(omega_matter-Cosmo.Omega_Matter) > 1e-3, 
			"Matter content of box is inconsistent with Omega_Matter in "
			"parameter file ! %g != %g ", omega_matter, Cosmo.Omega_Matter );

#endif // COMOVING

	return ;
}
