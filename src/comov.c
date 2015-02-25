#include "globals.h"
#include "kick.h"
#include "drift.h"

#ifdef COMOVING

#define TABLESIZE 1000

static float drift_factor_table[TABLESIZE] = { 0 };
static float kick_factor_table[TABLESIZE] = { 0 };

double Particle_Kick_Step(const int ipart)
{
	return 1;
}

double Particle_Drift_Step()
{


	return 1;
}



void Setup_Comoving()
{

	return ;
}




static float symplectic_comoving_drift_factor() // Quinn+ 1996
{
	return 1;
}

static float symplectic_comoving_kick_factor() // Quinn+ 1996
{
	return 1;
}



#endif // COMOVING
