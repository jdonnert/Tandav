#include "globals.h"
#include "timestep.h"
#include "drift.h"
#include "Gravity/gravity.h"

#ifndef COMOVING 
static double Particle_Drift_Step(const int ipart, const double time_next);
#endif

/* 
 * This is the drift part of the KDK scheme (Dehnen & Read 2012, Springel 05). 
 * As a snapshot time may not fall onto an integertime, we have to 
 * drift to the snapshot time, write the snapshot and then drift the 
 * remaining time to the next integertime. We use signals ... 
 */

void Drift_To_Sync_Point()
{
	Profile("Drift");

	#pragma omp for
	for (int i = 0; i < NParticle_Vectors; i++) {

		Float dt = Particle_Drift_Step(V.First[i], Time.Next);

		#pragma IVDEP
		for (int ipart = V.First[i]; ipart < V.Last[i]; ipart++) {

			P.Pos[0][ipart] += dt * P.Vel[0][ipart];
			P.Pos[1][ipart] += dt * P.Vel[1][ipart];
			P.Pos[2][ipart] += dt * P.Vel[2][ipart];
		}
	}

	if (!Sig.Domain_Update)
		Gravity_Tree_Update_Drift(Time.Step);

	Periodic_Constrain_Particles_To_Box();

	#pragma omp single
	{

	Int_Time.Current += Int_Time.Step;

	Time.Current = Integer_Time2Integration_Time(Int_Time.Current);

	Time.Step_Counter++;

	} // omp single

	Set_Current_Cosmology(Time.Current); // update immediately

	Profile("Drift");

	return;
}

/* 
 * Drift the system forward only to the snaptime the system is then 
 * NOT synchronized with the integer timeline. This is corrected during 
 * the next drift.
 */

void Drift_To_Snaptime()
{
	rprintf("\nDrift to next Shapshot Time %g -> %g \n", Time.Current, 
			Time.Next_Snap);

	const intime_t int_time_pos = Integration_Time2Integer_Time(Time.Next_Snap);
	
	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		Float dt = Particle_Drift_Step(ipart, Time.Next_Snap);

		P.Pos[0][ipart] +=	dt * P.Vel[0][ipart];
		P.Pos[1][ipart] +=	dt * P.Vel[1][ipart];
		P.Pos[2][ipart] +=	dt * P.Vel[2][ipart];

		P.Int_Time_Pos[ipart] = int_time_pos;
	}

	Periodic_Constrain_Particles_To_Box();

	#pragma omp single
	Time.Current = Time.Next_Snap;
	
	Set_Current_Cosmology(Time.Current); 

	return ;
}

/*
 * Return the amount of time since the last drift of particle ipart.
 */

#ifndef COMOVING 

static double Particle_Drift_Step(const int ipart, const double time_next)
{
	double time_part = Integer_Time2Integration_Time(P.Int_Time_Pos[ipart]);

	return time_next - time_part;
}

#endif // COMOVING


