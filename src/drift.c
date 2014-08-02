#include "globals.h"
#include "timestep.h"

/* This is the drift part of the KDK scheme (Dehnen & Read 2012, Springel 05). 
 * As a snapshot time may not fall onto an integertime, we have to 
 * drift to the snapshot time, write the snapshot and then drift the 
 * remaining time to the next integertime */

void Drift_To_Sync_Point() 
{
#ifdef COMOVING
	const float driftfac = Cosmo_Drift_Factor(Time.Current);
#else
	const float driftfac = 1;
#endif

	double dt = Time.Step;

	if (Sig.Synchronize_Drift) { // drift the rest to next integer time 
	
		dt = Integer2PhysicalTime(Time.IntNext) - Time.Current;

		Sig.Synchronize_Drift = false;
	}

	#pragma omp parallel for 
	for (int i = 0; i < NActiveParticles; i++) {

		int ipart = ActiveParticleList[i];
		if (ipart == 0)
			continue;

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0] * driftfac;
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1] * driftfac;
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2] * driftfac;
	}
	
	Time.IntCurrent += Time.IntStep;

	Time.Current = Integer2PhysicalTime(Time.IntCurrent);

	return;
}

/* drift the system forward only to the snaptime 
 * the system is then not synchronized with the 
 * integer timeline */

void Drift_To_Snaptime()
{
	const double dt = Time.NextSnap - Time.Current; // only drift this far

#ifdef COMOVING
	const float driftfac = Cosmo_Drift_Factor(Time.Current);
#else
	const float driftfac = 1;
#endif

	#pragma omp parallel for 
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {

		if (ipart == 0)
			continue;

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0] * driftfac;
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1] * driftfac;
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2] * driftfac;
	}

	Time.Current += dt;
	
	Time.NextSnap += Time.BetSnap;

	Sig.Synchronize_Drift = true; // signal Drift() to do the rest only

	return ;
}
