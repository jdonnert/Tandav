#include "unit.h"

const struct Code_Units Unit = {
	LENGTH2CGS,
	MASS2CGS,
	VELOCITY2CGS,
	LENGTH2CGS/VELOCITY2CGS, // time
	MASS2CGS*p2(VELOCITY2CGS), // energy
	MASS2CGS/p3(LENGTH2CGS) // density
};

void Init_Units()
{
	printf("Units: \n"
			"  Length   = %g cm\n"
			"  Mass     = %g g\n"
			"  Velocity = %g cm/s\n"
			"  Time     = %g s \n"
			"  Energy   = %g erg\n"
			"  Density  = %g g/cm^3 \n\n",
			Unit.Length, Unit.Mass, Unit.Velocity, Unit.Time, Unit.Energy, 
			Unit.Density);

	return ;
}

/* Mostly for reference purposes, these functions convert internal units 
 * to physical cgsm units. In case of comoving coordinates, they also 
 * remove the a factors. */

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

double Acceleration_Physical(const int ipart)
{
	double grav_accel = sqrt( p2(P.Grav_Acc[0][ipart]) 
			+ p2(P.Grav_Acc[1][ipart]) + p2(P.Grav_Acc[2][ipart]) );

#ifdef COMOVING
	grav_accel *= Cosmo.Grav_Accel_Factor;
#endif

	double hydro_accel = 0;

#ifdef COMOVING
	hydro_accel *= Cosmo.Hydro_Accel_Factor;
#endif

	return grav_accel + hydro_accel;
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

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
