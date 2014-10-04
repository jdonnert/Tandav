#include "globals.h"
#include "timestep.h"

/* 
 * This is the Kick part of the KDK scheme. We update velocities from 
 * accelerations, but kick only for half a timestep 
 */

void Kick_First_Halfstep() 
{
	Profile("First Kick");

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

#ifdef COMOVING 
		double dt = 0.5 * Cosmo_Kick_Factor(Sim.Current_Time); // Quinn+97
#else
		double dt = 0.5 * (Time.Next - time_part);
#endif // COMOVING

		P[ipart].Vel[0] += dt * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt * P[ipart].Acc[1]; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2];
	}
	
	Profile("First Kick");

	return ;
}

void Kick_Second_Halfstep() 
{
	Profile("Second Kick");

	#pragma omp  for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

#ifdef COMOVING 
		double dt = 0.5 * Cosmo_Kick_Factor(Sim.Current_Time); // Quinn+97
#else
		double dt = 0.5 * (Time.Next - time_part);
#endif // COMOVING

		P[ipart].Vel[0] += dt * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt * P[ipart].Acc[1]; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2];

		P[ipart].Int_Time_Pos = Int_Time.Next;
	}
	
	Profile("Second Kick");

	return ;
}
