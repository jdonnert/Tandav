#include "globals.h"
#include "timestep.h"

/* 
 * This is the Kick part of the KDK scheme. We update velocities from forces,
 * but kick only for half a timestep 
 */

void Kick_Halfstep() 
{
#ifdef COMOVING 
	const float kick_fac = Cosmo_Kick_Factor(Sim.CurrentTime); // Quinn+97
#else
	const float kick_fac = 1;
#endif // COMOVING
	
	#pragma omp parallel for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		float dt = 0.5 * Timebin2Timestep(P[ipart].TimeBin);

		P[ipart].Vel[0] += dt * P[ipart].Acc[0] * kick_fac;
		P[ipart].Vel[1] += dt * P[ipart].Acc[1] * kick_fac; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2] * kick_fac;
	}
	
	return ;
}

