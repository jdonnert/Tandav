#include "globals.h"
#include "timestep.h"

/* 
 * This is the Kick part of the KDK scheme. We update velocities from forces,
 * but kick only for half a timestep 
 */

void Kick_Halfstep() 
{
#ifdef COMOVING 
	const float kick_factor = Cosmo_Kick_Factor(Sim.Current_Time); // Quinn+97
#else
	const float kick_factor = 1;
#endif // COMOVING
	
	#pragma omp parallel for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		float dt = 0.5 * Timebin2Timestep(P[ipart].Time_Bin);

		P[ipart].Vel[0] += dt * P[ipart].Acc[0] * kick_factor;
		P[ipart].Vel[1] += dt * P[ipart].Acc[1] * kick_factor; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2] * kick_factor;
	}
	
	return ;
}

