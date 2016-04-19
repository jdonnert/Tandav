#include "globals.h"
#include "timestep.h"
#include "drift.h"
#include "Gravity/gravity.h"

/* 
 * This is the drift part of the KDK scheme (Dehnen & Read 2012, Springel 05). 
 * As a snapshot time may not fall onto an integertime, we have to 
 * drift to the snapshot time, write the snapshot and then drift the 
 * remaining time to the next integertime. We use signals ... 
 */

void Drift_To_Sync_Point()
{
	Profile("Drift");

	#pragma omp single //for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		intime_t it_step = Timebin2It_Timestep(P.Time_Bin[ipart]);

		intime_t it_curr = P.It_Drift_Pos[ipart];
		intime_t it_next = it_curr + it_step;

		Assert(it_next <= Int_Time.End, 
				"overstepped ipart=%d, curr=%u next=%u max=%u"
				"IT.curr=%u IT.next=%u ", 
				ipart, it_curr, it_next, Int_Time.End, Int_Time.Current, 
				Int_Time.Next);

		double dt = Particle_Drift_Step(ipart, it_curr, it_next);

		P.Pos[0][ipart] += dt * P.Vel[0][ipart];
		P.Pos[1][ipart] += dt * P.Vel[1][ipart];
		P.Pos[2][ipart] += dt * P.Vel[2][ipart];

		P.It_Drift_Pos[ipart] += it_step;
	}

	if (!Sig.Domain_Update)
		Gravity_Tree_Update_Drift(Time.Step);

	Periodic_Constrain_Particles_To_Box(); // PERIODIC

	#pragma omp single
	{

	Int_Time.Current += Int_Time.Step;
	Int_Time.Next += Int_Time.Step;

	Time.Current = Integer_Time2Integration_Time(Int_Time.Current);
	Time.Next = Integer_Time2Integration_Time(Int_Time.Next);

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

	const intime_t it_next = Integration_Time2Integer_Time(Time.Next_Snap);
	
	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		intime_t it_curr = P.It_Drift_Pos[ipart];

		double dt = Particle_Drift_Step(ipart, it_curr, it_next);

		P.Pos[0][ipart] +=	dt * P.Vel[0][ipart];
		P.Pos[1][ipart] +=	dt * P.Vel[1][ipart];
		P.Pos[2][ipart] +=	dt * P.Vel[2][ipart];

		P.It_Drift_Pos[ipart] = it_next;
		P.Time_Bin[ipart] = log2(it_next - it_curr); // dt_it -> timebin
	}

	Periodic_Constrain_Particles_To_Box();

	#pragma omp single
	Time.Current = Time.Next_Snap;
	
	Set_Current_Cosmology(Time.Current); 

	return ;
}

/*
 * Return the amount of real time between two points on the integer timeline.
 * This is a 2^X of the smallest increment representable on the integer time-
 * line, when we are not in comoving coordinates. 
 */

#ifndef COMOVING 

double Particle_Drift_Step(const int ipart ,
		const intime_t it_curr, const intime_t it_next)
{
	double t_curr = Integer_Time2Integration_Time(it_curr);
	double t_next = Integer_Time2Integration_Time(it_next);

	return t_next - t_curr;
}

#endif // COMOVING


