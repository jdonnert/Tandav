#include "globals.h"
#include "timestep.h"

/* 
 * This is the Kick part of the KDK scheme. We update velocities from 
 * accelerations, but kick only for half a timestep 
 */

void Kick_Halfstep() 
{
	#pragma omp parallel for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		float dt_native = 0.5 * Timebin2Timestep(P[ipart].Time_Bin);

#ifdef COMOVING 
		float dt_grav = 0.5 * Cosmo_Kick_Factor(Sim.Current_Time); // Quinn+97
#else
		float dt_grav = dt_native;
#endif // COMOVING

		P[ipart].Vel[0] += dt_grav * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt_grav * P[ipart].Acc[1]; 
		P[ipart].Vel[2] += dt_grav * P[ipart].Acc[2];
	}
	
	return ;
}

