#include "globals.h"
#include "timestep.h"

struct CurrentCosmologyInCodeUnits Cosmo = { 
	HUBBLE_CONST, OMEGA_LAMBDA, 
	(OMEGA_0 - OMEGA_LAMBDA - OMEGA_RAD)*p2(HUBBLE_CONST/100), // Omega_M
	OMEGA_BARYON * p2(HUBBLE_CONST/100), OMEGA_0, 
	OMEGA_RAD / p2(HUBBLE_CONST/100), 
	3/8/PI/(GRAVITATIONAL_CONST/p3(VELOCITY2CGS)/TIME2CGS*MASS2CGS)
		*p2(HUBBLE_CONST), 	// rho_crit0, can't use Const.Gravity here :-(
	0 };

/*
 * This updates the Cosmo structure to the current Time.Current == a
 */

void Set_Current_Cosmology()
{
#ifndef COMOVING
	Assert(false, "Recompile with COMOVING defined");
#endif

	const double a = Time.Current; // yes this will be everywhere in the code
	
	Cosmo.Expansion_Factor = a;

	Cosmo.Redshift = 1/a - 1;
	Cosmo.Hubble_Parameter = Hubble_Parameter(a);
	Cosmo.Critical_Density = Critical_Density(a);

	return ;
}

/* 
 * These functions implements the cosmological background evolution 
 * see also Peebles 1980
 */

double Hubble_Parameter(const double a) // H(a) Mo, v.d.Bosch & White 3.74
{
	return Cosmo.Hubble_Constant * E_Hubble(a); 
}

double E_Hubble(const double a) // E(a) Mo, v.d.Bosch & White 3.75
{
	return sqrt(OMEGA_LAMBDA + (1-OMEGA_0)/(a*a) 
			+ Cosmo.Omega_Matter/(a*a*a) + OMEGA_RAD/(a*a*a*a));
}

double Critical_Density(double hubble_parameter) // Mo, v.d.Bosch & White 3.63
{
	return 3.0/8.0/PI/Const.Gravity * p2(hubble_parameter);
}
