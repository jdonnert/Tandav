#include "globals.h"
#include "timestep.h"
#include "kick.h"
#include "Gravity/gravity.h"

#ifndef COMOVING 
double Particle_Kick_Step(const int ipart, const double time_next);
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

		double dt = 0.5 * Particle_Kick_Step(ipart, Time.Next);

		P.Vel[0][ipart] += dt * P.Acc[0][ipart];
		P.Vel[1][ipart] += dt * P.Acc[1][ipart];
		P.Vel[2][ipart] += dt * P.Acc[2][ipart];

		if (!Sig.Domain_Update)
			Gravity_Tree_Update_Kicks(ipart, dt); // GRAVITY_TREE
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

		double dt = 0.5 * Particle_Kick_Step(ipart, Time.Next);

		P.Vel[0][ipart] += dt * P.Acc[0][ipart];
		P.Vel[1][ipart] += dt * P.Acc[1][ipart];
		P.Vel[2][ipart] += dt * P.Acc[2][ipart];

		P.Int_Time_Pos[ipart] = Int_Time.Next;

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

double Particle_Kick_Step(const int ipart, const double time_next)
{
	const intime_t intime_pos = P.Int_Time_Pos[ipart];

	double time_part = Integer_Time2Integration_Time(intime_pos);

	return time_next - time_part;
}

#endif // ! COMOVING
