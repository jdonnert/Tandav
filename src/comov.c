#ifdef COMOVING
#include "globals.h"

#define TABLESIZE 1000

static float drift_factor_table[TABLESIZE] = { 0 };
static float kick_factor_table[TABLESIZE] = { 0 };

static float symplectic_comoving_drift_factor() // Quinn+ 1996
{
	return 1;
}

static float symplectic_comoving_kick_factor() // Quinn+ 1996
{
	return 1;
}

#endif // COMOVING
