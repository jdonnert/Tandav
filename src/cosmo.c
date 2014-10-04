#include "globals.h"

struct Cosmology_Infos Cosmo = { 0 };
#pragma omp threadprivate(Cosmo)

double Cosmo_Drift_Factor(double a)
{
	return 1;
}

double Cosmo_Kick_Factor(double a)
{
	return 1;
}

double Hubble_Function(double a) // i.e. Mo, van den Bosch, White 2010
{
	return  1;
}

double Critical_Density(double a)
{
	return 1;
}

void Init_Cosmology()
{

	return ;
}

