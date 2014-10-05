#include "globals.h"
#include "timestep.h"

/* 
 * This is the drift part of the KDK scheme (Dehnen & Read 2012, Springel 05). 
 * As a snapshot time may not fall onto an integertime, we have to 
 * drift to the snapshot time, write the snapshot and then drift the 
 * remaining time to the next integertime 
 */

void Drift_To_Sync_Point() 
{
	Profile("Drift");

	rprintf("Drift to next Sync Point ... ");
 
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

		double dt = Time.Next - time_part;
		
		P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];

	}
	
	#pragma omp single 
	{

	Int_Time.Current += Int_Time.Step;

	Time.Current = Integer2Physical_Time(Int_Time.Current);

	mprintf("done \n");
	
	}

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
	rprintf("Drift to next Shapshot Time %g \n", Time.Next_Snap);

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		double dt = Time.Next_Snap 
			- Integer2Physical_Time(P[ipart].Int_Time_Pos);

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];
	}
	
	#pragma omp single
	Time.Current = Time.Next_Snap;

	return ;
}
