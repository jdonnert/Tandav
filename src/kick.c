#include "kick.h"

/*
 * This is the Kick part of the KDK scheme. We update velocities from
 * accelerations, but kick only for half a timebin. If we use the tree, the
 * nodes are kicked as well.
 */

void Kick_First_Halfstep()
{
	Profile("First Kick");

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		intime_t it_step = Timebin2It_Timestep(P.Time_Bin[ipart]);

		intime_t it_curr = P.It_Kick_Pos[ipart];
		intime_t it_next = it_curr + (it_step >> 1);

		Float dt = Particle_Kick_Step(it_curr, it_next);

		P.Vel[0][ipart] += dt * P.Acc[0][ipart];
		P.Vel[1][ipart] += dt * P.Acc[1][ipart];
		P.Vel[2][ipart] += dt * P.Acc[2][ipart];

		P.It_Kick_Pos[ipart] += (it_step >> 1);
	} // for i

	Gravity_Tree_Update_Kicks(); // GRAVITY_TREE

	Profile("First Kick");

	return ;
}

void Kick_Second_Halfstep()
{
	Profile("Second Kick");

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		intime_t it_step = Timebin2It_Timestep(P.Time_Bin[ipart]);

		intime_t it_curr = P.It_Kick_Pos[ipart];
		intime_t it_next = P.It_Kick_Pos[ipart] + (it_step >> 1);

		Float dt = Particle_Kick_Step(it_curr, it_next);

		P.Vel[0][ipart] += dt * P.Acc[0][ipart];
		P.Vel[1][ipart] += dt * P.Acc[1][ipart];
		P.Vel[2][ipart] += dt * P.Acc[2][ipart];

		P.It_Kick_Pos[ipart] += (it_step >> 1);
	}

	Gravity_Tree_Update_Kicks(); // GRAVITY_TREE

	Profile("Second Kick");

	return ;
}

/*
 * Return the amount of real time between two points on the integer timeline.
 */

#ifndef COMOVING

double Particle_Kick_Step(const intime_t it_curr, const intime_t it_next)
{
	double t_curr = Integer_Time2Integration_Time(it_curr);
	double t_next = Integer_Time2Integration_Time(it_next);

	return t_next - t_curr;
}

#endif // ! COMOVING
