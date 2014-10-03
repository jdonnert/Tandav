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
 
	const double t_next = Integer2Physical_Time(Int_Time.Next);
	const double dt_max = t_next - Time.Current;

	#pragma omp parallel for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		//double dt = t_next - Integer2Physical_Time(P[ipart].Int_Time_Pos);
	double dt = Timebin2Timestep(P[ipart].Time_Bin);

		//dt = fmin(dt, dt_max);

		//dt = fmin(dt, t_next - Time.Current);
//	if (dt != dt2)
//		printf("i=%d dt=%g dt2=%g cur=%g nex=%g it=%g\n", ipart, dt, dt2, Time.Current, t_next,
//				Integer2Physical_Time(P[ipart].Int_Time_Pos));

		P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];

		P[ipart].Int_Time_Pos = Int_Time.Next;
	}
	
	Int_Time.Current += Int_Time.Step;

	Time.Current = Integer2Physical_Time(Int_Time.Current);

	rprintf("done \n");

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
{return ;
	rprintf("Drift to next Shapshot Time ...");

	#pragma omp parallel for
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];

		double dt = Time.Next_Snap 
			- Integer2Physical_Time(P[ipart].Int_Time_Pos);

	 	P[ipart].Pos[0] += 	dt * P[ipart].Vel[0];
		P[ipart].Pos[1] += 	dt * P[ipart].Vel[1];
		P[ipart].Pos[2] += 	dt * P[ipart].Vel[2];
	}

	Time.Current = Time.Next_Snap;
	
	rprintf("done \n");

	return ;
}
