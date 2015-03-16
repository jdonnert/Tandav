#include "globals.h"
#include "timestep.h"
#include "kick.h"
#include "Gravity/gravity.h"

#ifndef COMOVING 
double Particle_Kick_Step(const int ipart);
#endif

/* 
 * This is the Kick part of the KDK scheme. We update velocities from 
 * accelerations, but kick only for half a timestep. If we use the tree, the
 * nodes are kicked as well.
 */

void Kick_First_Halfstep()
{
	Profile("First Kick");

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		double dt = Particle_Kick_Step(ipart);

		P[ipart].Vel[0] += dt * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt * P[ipart].Acc[1];
		P[ipart].Vel[2] += dt * P[ipart].Acc[2];

		if (!Sig.Domain_Update)
			Gravity_Tree_Update_Kicks(ipart, dt);
	}

#pragma omp barrier

	Profile("First Kick");

	return ;
}

void Kick_Second_Halfstep()
{
	Profile("Second Kick");

	#pragma omp  for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		double dt =  0.5 * Particle_Kick_Step(ipart);

		P[ipart].Vel[0] += dt * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt * P[ipart].Acc[1];
		P[ipart].Vel[2] += dt * P[ipart].Acc[2];

		P[ipart].Int_Time_Pos = Int_Time.Next;

		if (!Sig.Domain_Update)
			Gravity_Tree_Update_Kicks(ipart, dt); 
	}

	Profile("Second Kick");

	return ;
}

/*
 * Return half the amount of time since the last kick of particle ipart. For
 * cosmological time integration this is the kick factor in comov.c
 */

#ifndef COMOVING 
double Particle_Kick_Step(const int ipart)
{
	double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

	return (Time.Next - time_part);
}
#endif
