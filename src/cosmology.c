#include "globals.h"
#include "timestep.h"

struct Current_Cosmology_In_Code_Units Cosmo = {
	HUBBLE_CONST * 1e5 / (1e3*KPC2CGS) * (LENGTH2CGS/VELOCITY2CGS), // H0
	OMEGA_LAMBDA,
	(OMEGA_0 - OMEGA_LAMBDA - OMEGA_RAD)*p2(HUBBLE_CONST/100.0), // Omega_M
	OMEGA_BARYON * p2(HUBBLE_CONST/100.0),
	OMEGA_0,
	OMEGA_RAD / p2(HUBBLE_CONST/100.0),
	3.0/8.0/PI / (GRAVITATIONAL_CONST/p3(VELOCITY2CGS)
			/(LENGTH2CGS/VELOCITY2CGS)*MASS2CGS)
			*p2(HUBBLE_CONST*1e5/(1e3*KPC2CGS)), // rho0_crit	
	0 // the rest is set in "Set_Current_Cosmology()"
};

#ifdef COMOVING

/*
 * This updates the Cosmo structure to the current Time.Current == a
 */

void Set_Current_Cosmology()
{
	const double a = Time.Current; // yes this will be everywhere in the code

	Cosmo.Expansion_Factor = a; // just to be clear

	Cosmo.Redshift = 1/a - 1;
	Cosmo.Hubble_Parameter = Hubble_Parameter(a);
	Cosmo.Critical_Density = Critical_Density(a);

	return ;
}

/* 
 * These functions implement the cosmological background evolution in
 * code units, see also Peebles 1980, Mo, v.d.Bosch & White Eq. 3.74/5.
 */

double Hubble_Parameter(const double a) // H(a) = H0 * E(a), (Eq 3.74)
{
	return Cosmo.Hubble_Constant * E_Hubble(a);
}

double E_Hubble(const double a) // E(a), (Eq 3.75)
{
	return sqrt(OMEGA_LAMBDA + (1-OMEGA_0)/(a*a)
			+ Cosmo.Omega_Matter/(a*a*a) + OMEGA_RAD/(a*a*a*a));
}

double Critical_Density(double hubble_parameter) // Mo, v.d.Bosch & White 3.63
{
	return 3.0/8.0/PI/Const.Gravity * p2(hubble_parameter);
}

void Setup_Cosmology()
{
	const double h0_cgs = HUBBLE_CONST * 1e5 / (1e3*KPC2CGS);

	rprintf("Cosmological Model: \n"
			"   h_0          = %4g, Omega_0      = %4g\n"
			"   Omega_Lambda = %4g, Omega_Matter = %4g\n"
			"   Omega_Baryon = %4g, Omega_Rad    = %4g\n"
			"   rho_crit_0   = %4g g/cm^3\n", 
			HUBBLE_CONST/100, Cosmo.Omega_0, Cosmo.Omega_Lambda, 
			Cosmo.Omega_Matter, Cosmo.Omega_Baryon, Cosmo.Omega_Rad, 
			3.0/8.0/PI/GRAVITATIONAL_CONST*p2(h0_cgs));
printf("%g \n", Cosmo.Hubble_Constant);
	return ;
}


#endif // COMOVING
