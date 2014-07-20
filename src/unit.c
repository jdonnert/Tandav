/* Conversion from code to CGSM or code units */

#include "globals.h"

const struct Units Unit = { 	
	LENGTH2CGS, 
	MASS2CGS, 
	VELOCITY2CGS, 
	LENGTH2CGS/VELOCITY2CGS, // time
	MASS2CGS*p2(VELOCITY2CGS) // energy
};

/* Get derived quantities in code units */
double Pressure(const int ipart)
{
	return 1; //Gas[ipart].Entropy * pow(Gas[ipart].Rho, Const.Adiabatic_Index);
}

double Internal_Energy(const int ipart)
{
	float density = 1; //Gas[ipart].Density;

#ifdef COMOVOING
	density /= p3(Time.Current);
#endif
	
	double u = 1 // Gas[ipart].Entropy / (Const.Adiabatic_Index-1) 
		* pow(density, (Const.Adiabatic_Index-1));

	return u;
}

double Temperature(const int ipart)
{
	return 1;
}

/* These functions convert internal units to physical cgsm units */
double Position_Cgs(const float x)
{
	double x_cgs = x * Unit.Length;

#ifdef COMOVING
	x_cgs *= Current.Time / Cosmo.HubbleParam; 
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
	mass_cgs /= Cosmo.HubbleParam; 
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
		/ p2(Cosmo.HubbleParam);
#endif

	return press_cgs;
}

double Thermal_Energy_Density_Cgs(const int ipart)
{
	return 1;
}
