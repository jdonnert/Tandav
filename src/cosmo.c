#include "globals.h"
#include "proto.h"


static float Hubble_Function () // i.e. Mo, van den Bosch, White 2010
{

	return ;
}

#define TABLESIZE 1000

static float drift_factor_table[TABLESIZE] = { 0 };
static float kick_factor_table[TABLESIZE] = { 0 };

static float symplectic_drift_factor() // Quinn+ 1996
{
	return 1;
}

static float symplectic_kick_factor() // Quinn+ 1996
{
	return 1;
}

void Fill_Comoving_Factor_Tables()
{

	return ;
}

#undef TABLESIZE
