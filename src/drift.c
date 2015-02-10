#include "globals.h"
#include "timestep.h"
#include "drift.h"
#include "Gravity/gravity.h"

/* 
 * This is the drift part of the KDK scheme 
 * (Dehnen & Read 2012, Springel 05). 
 * As a snapshot time may not fall onto an integertime, we have to 
 * drift to the snapshot time, write the snapshot and then drift the 
 * remaining time to the next integertime. We use signal 
 */

static double time_snap = 0; 

void Drift_To_Sync_Point() 
{
	Profile("Drift");
	
	if (Sig.Drifted_To_Snaptime) { // handle out of sync integer timeline

		Sig.Drifted_To_Snaptime = false;
	
		#pragma omp single
		time_snap = Time.Current;
	}

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

		time_part = fmax(time_part, time_snap);

		double dt = Time.Next - time_part;
		
		P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];

	}

	Constrain_Particles_To_Box(); // PERIODIC

#ifdef GRAVITY_TREE
	if (! Sig.Domain_Update)
		Gravity_Tree_Update_Drift(Time.Step); 
#endif

	#pragma omp single 
	{

	Int_Time.Current += Int_Time.Step;

	Time.Current = Integer2Physical_Time(Int_Time.Current);

	} // omp single

	#pragma omp barrier

	Profile("Drift");
	
	return;
}

/* 
 * drift the system forward only to the snaptime 
 * the system is then NOT synchronized with the 
 * integer timeline. This is corrected during the next
 * drift.
 */

void Drift_To_Snaptime()
{
	rprintf("\nDrift to next Shapshot Time %g \n", Time.Next_Snap);

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		double dt = Time.Next_Snap 
			- Integer2Physical_Time(P[ipart].Int_Time_Pos);

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];
	}

	Constrain_Particles_To_Box(); // PERIODIC

	Sig.Drifted_To_Snaptime = true;

	#pragma omp single
	Time.Current = Time.Next_Snap;

	return ;
}

#ifdef PERIODIC
static void Constrain_Particles_To_Box()
{
	const Float boxsize[3] = {Sim.Boxsize[0], Sim.Boxsize[1], Sim.Boxsize[2]};

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		while (P[ipart].Pos[0] < 0)
			P[ipart].Pos[0] += boxsize[0];
		
		while (P[ipart].Pos[0] >= boxsize[0])
			P[ipart].Pos[0] -= boxsize[0];

		while (P[ipart].Pos[1] < 0)
			P[ipart].Pos[1] += boxsize[1];
		
		while (P[ipart].Pos[1] >= boxsize[1])
			P[ipart].Pos[1] -= boxsize[1];

		while (P[ipart].Pos[2] < 0)
			P[ipart].Pos[2] += boxsize[2];
		
		while (P[ipart].Pos[2] >= boxsize[2])
			P[ipart].Pos[2] -= boxsize[2];
	} // for i

	return ;
}
#endif
