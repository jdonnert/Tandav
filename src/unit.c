/* 
 * Conversion from code to CGSM or code units 
 */

#include "globals.h"

const struct Units Unit = { 	
	LENGTH2CGS, 
	MASS2CGS, 
	VELOCITY2CGS, 
	LENGTH2CGS/VELOCITY2CGS, // time
	MASS2CGS*p2(VELOCITY2CGS) // energy
};

/* 
 * These functions convert internal units to physical cgsm units 
 */

double Position_Cgs(const float x)
{
	double x_cgs = x * Unit.Length;

#ifdef COMOVING
	x_cgs *= Current.Time / Cosmo.Hubble_Param; 
#endif

	return x_cgs;
}

double Velocity_Cgs(const float v)
{
	double v_cgs = v * Unit.Velocity;

#ifdef COMOVING
	v_cgs *= Time.Current;
#endif
	
	return v_cgs;
}

double Mass_Cgs(const float mass)
{
	double mass_cgs = mass * Unit.Mass;

#ifdef COMOVING
	mass_cgs /= Cosmo.Hubble_Param; 
#endif

	return mass_cgs;
}	

double Density_Cgs(const float rho)
{
	double rho_cgs = rho * Unit.Mass / p3(Unit.Length);

#ifdef COMOVING
	rho_cgs *= p2(Cosmo.HubbleParam) / p3(Time.Current);
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
	press_cgs *= pow(Time.Current, -3*Const.Adiabatic_Index) 
		/ p2(Cosmo.Hubble_Param);
#endif

	return press_cgs;
}

double Thermal_Energy_Density_Cgs(const int ipart)
{
	return 1;
}
