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

	Setup_Domain_Decomposition();
	
	Setup_Gravity_Tree(); // GRAVITY_TREE

	/* Add yours above */
	
	Compute_Global_Simulation_Properties();

	sanity_check_simulation_setup();
	
	Print_Memory_Usage();

	Profile("Setup");

	return ;
}


static void sanity_check_simulation_setup()
{

#ifdef COMOVING

	double rho_crit = 3 * p2(Cosmo.Hubble_Constant)/(8 * PI * Const.Gravity);
	double omega_matter = Sim.Total_Mass / p3(Sim.Boxsize[0]) / rho_crit;

	Warn(fabs(omega_matter-Cosmo.Omega_Matter) > 1e-3, 
			"Matter content of box is inconsistent with Omega_Matter in "
			"parameter file. \n O_M ICs %g : Code %g. Did you set HUBBLE_CONST"
			" = %g correctly ?", 
			omega_matter, Cosmo.Omega_Matter, HUBBLE_CONST );

#endif // COMOVING

#ifdef GRAVITY

	double volume = Sim.Boxsize[0]*Sim.Boxsize[1]*Sim.Boxsize[2]; 
	double mean_part_dist = pow(volume / Sim.Npart_Total, 1.0/3.0);
	double eps_std = mean_part_dist/7.0;

	Warn(Param.Grav_Softening[1] > eps_std*10, "GravSoftening seems large, "
			"have %g recommend %g", Param.Grav_Softening[1], eps_std);

	Warn(Param.Grav_Softening[1] < eps_std/10, "GravSoftening seems small, "
			"have %g recommend %g", Param.Grav_Softening[1], eps_std);

#endif // GRAVITY

	return ;
}
