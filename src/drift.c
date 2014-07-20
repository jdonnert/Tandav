#include "globals.h"
#include "timestep.h"

void Drift() 
{
#ifdef COMOVING
	const float driftfac = Cosmo_Drift_Factor(Time.Current);
#else
	const float driftfac = 1;
#endif

	const float dt = Time.Step;

#pragma omp parallel for 
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {

		if (ipart == 0)
			continue;

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0] * driftfac;
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1] * driftfac;
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2] * driftfac;
	}
	
	Time.Current += dt; // we have advanced the system !

	return;
}

