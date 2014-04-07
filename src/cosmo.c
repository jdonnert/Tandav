#include "globals.h"

struct Cosmology_Infos Cosmo = { 0 };

double Cosmo_Drift_Factor(float a)
{
	return 1;
}
double Cosmo_Kick_Factor(float a)
{
	return 1;
}
float Hubble_Function (float a) // i.e. Mo, van den Bosch, White 2010
{
	return  1;
}

float Critical_Density (float a)
{
	return 1;
}

void Init_Cosmology ()
{

	return ;
}

#undef TABLESIZE
