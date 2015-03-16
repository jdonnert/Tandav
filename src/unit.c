/* 
 * Conversion from code to CGSM units 
 */

#include "globals.h"

const struct Code_Units Unit = { 	
	LENGTH2CGS, 
	MASS2CGS, 
	VELOCITY2CGS, 
	LENGTH2CGS/VELOCITY2CGS, // time
	MASS2CGS*p2(VELOCITY2CGS) // energy
};

/* 
 * These functions convert internal units to physical cgsm units. In case of 
 * comoving coordinates, they also remove the a factors.
 */

double Position_Cgs(const float x)
{
	double x_cgs = x * Unit.Length;

#ifdef COMOVING
	x_cgs *= Cosmo.Expansion_Factor; 
#endif

	return x_cgs;
}

double Velocity_Cgs(const float v)
{
	double v_cgs = v * Unit.Velocity;

#ifdef COMOVING
	v_cgs *= p2(Cosmo.Expansion_Factor);
#endif
	
	return v_cgs;
}

double Mass_Cgs(const float mass)
{
	double mass_cgs = mass * Unit.Mass;

	return mass_cgs;
}	

double Density_Cgs(const float rho)
{
	double rho_cgs = rho * Unit.Mass / p3(Unit.Length);

#ifdef COMOVING
	rho_cgs /= p3(Cosmo.Expansion_Factor);
#endif

	return rho_cgs;
}

double Number_Density_Cgs(const float rho)
{ 
	return rho * ( (3*HYDROGEN_FRACTION + 1) / (4 * PROTON_MASS));
}

double Pressure_Cgs(const float press)
{
	double press_cgs = press * Unit.Energy / p3(Unit.Length);

#ifdef COMOVING
	press_cgs *= pow(Cosmo.Expansion_Factor, -3*Const.Adiabatic_Index);
#endif

	return press_cgs;
}

double Thermal_Energy_Density_Cgs(const int ipart)
{
	return 1;
}
