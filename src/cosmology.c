#include "globals.h"
#include "timestep.h"

#pragma omp threadprivate(Cosmo)
extern struct Current_Cosmology_In_Code_Units Cosmo = {
	HUBBLE_CONST * KM2CGS/MPC2CGS * LENGTH2CGS/VELOCITY2CGS, // H0
	OMEGA_LAMBDA,
	OMEGA_MATTER,
	OMEGA_BARYON * p2(HUBBLE_CONST/100.0),
	OMEGA_MATTER + OMEGA_LAMBDA + OMEGA_RAD, // Mo+, eq. 3.72
	OMEGA_RAD / p2(HUBBLE_CONST/100.0),
	3.0/8.0/PI / (GRAVITATIONAL_CONST/p3(VELOCITY2CGS) // rho0_crit
			/(LENGTH2CGS/VELOCITY2CGS)*MASS2CGS)
			*p2(HUBBLE_CONST* KM2CGS/MPC2CGS * LENGTH2CGS/VELOCITY2CGS), 
	0 // the rest is done in "Set_Current_Cosmology()"
};

#ifdef COMOVING

void Setup_Cosmology()
{
	const double h0_cgs = HUBBLE_CONST * KM2CGS / MPC2CGS;

	rprintf("Cosmological background model: \n"
			"   h_0          = %6.3g, Omega_0      = %6.3g \n"
			"   Omega_Lambda = %6.3g, Omega_Matter = %6.3g \n"
			"   Omega_Baryon = %6.3g, Omega_Rad    = %6.3g \n"
			"   rho_crit_0   = %6.3g g/cm^3\n",
			HUBBLE_CONST/100, Cosmo.Omega_0, Cosmo.Omega_Lambda,
			Cosmo.Omega_Matter, Cosmo.Omega_Baryon, Cosmo.Omega_Rad,
			3.0/8.0/PI/GRAVITATIONAL_CONST*p2(h0_cgs));

	return ;
}

/*
 * This updates the variable parts of the Cosmo structure to the current 
 * expansion factor. 
 */

void Set_Current_Cosmology()
{
	const double a = Time.Current;

	Cosmo.Expansion_Factor = a;
	Cosmo.Sqrt_Expansion_Factor = sqrt(a);

	Cosmo.Redshift = 1/a - 1;
	Cosmo.Hubble_Parameter = Hubble_Parameter(a);
	Cosmo.Critical_Density = Critical_Density(a);

	Cosmo.Grav_Accel_Factor = 1/p2(a);
	Cosmo.Hydro_Accel_Factor = 1/pow(a, 3*(ADIABATIC_INDEX_MONOATOMIC_GAS - 2));
	Cosmo.Press_Factor = pow(a, 3*(ADIABATIC_INDEX_MONOATOMIC_GAS - 1));

	return ;
}

/* 
 * These functions implement the cosmological background evolution in
 * code units, see also Peebles 1980, Mo, v.d.Bosch & White Eq. 3.74/5.
 */

double Hubble_Parameter(const double a) // H(a) = H0 * E(a), Mo+ eq 3.74
{
	return Cosmo.Hubble_Constant * E_Hubble(a);
}

double E_Hubble(const double a) // E(a), Mo+ eq 3.75
{
	return sqrt(Cosmo.Omega_Lambda + (1.0 - Cosmo.Omega_0)/(a*a)
			+ Cosmo.Omega_Matter/(a*a*a) + OMEGA_RAD/(a*a*a*a));
}

double Critical_Density(double hubble_parameter) // Mo+ eq. 3.63
{
	return 3.0/(8.0*PI*Const.Gravity) * p2(hubble_parameter);
}

#endif // COMOVING
