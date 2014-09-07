#include "globals.h"
#include "kick.h"
#include "force.h"
#include "timestep.h"

/* This is the Kick part of the KDK scheme. We update velocities from forces,
 * but kick only for half a timestep */

void Kick_Halfstep() 
{
#ifdef COMOVING 
	const float driftfac = Cosmo_Kick_Factor(Sim.CurrentTime); // Quinn+97
#else
	const float driftfac = 1;
#endif // COMOVING
	
	for (int i = 0; i < NActiveParticles; i++) {

		int ipart = ActiveParticleList[i];
		
		float dt = 0.5 * Timebin2Timestep(P[ipart].TimeBin);

		P[ipart].Vel[0] += dt * P[ipart].Acc[0] * driftfac;
		P[ipart].Vel[1] += dt * P[ipart].Acc[1] * driftfac; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2] * driftfac;
	}
	
	return ;
}

