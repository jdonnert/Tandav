/* Unit conversion from code to CGSM units */
#include "globals.h"

const struct Units_In_Cgs Unit = { 	
	LENGTH2CGS, 
	MASS2CGS, 
	VELOCITY2CGS, 
	LENGTH2CGS/VELOCITY2CGS, // time
	MASS2CGS*p2(VELOCITY2CGS) // energy
};

double Density_Cgs(const int ipart)
{
	return 1;
};

double Temperature_Cgs(const int ipart)
{
	return 1;
};
