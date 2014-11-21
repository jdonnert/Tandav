#include "globals.h"
#include "timestep.h"

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

		double time_part = Integer2Physical_Time(P[ipart].Int_Time_Pos);

#ifdef COMOVING 
		double dt = 0.5 * Cosmo_Kick_Factor(Sim.Current_Time); // Quinn+97
#else
		double dt = 0.5 * (Time.Next - time_part);
#endif // COMOVING

		P[ipart].Vel[0] += dt * P[ipart].Acc[0];
		P[ipart].Vel[1] += dt * P[ipart].Acc[1]; 
		P[ipart].Vel[2] += dt * P[ipart].Acc[2];

#ifdef GRAVITY_TREE
		Float dv[3] = { dt*P[ipart].Acc[0], dt*P[ipart].Acc[1], 
			dt*P[ipart].Acc[2] };

		Gravity_Tree_Update_Kicks(dv, P[ipart].Tree_Parent); // kick tree nodes
#endif // GRAVITY_TREE
	}
	
	Gravity_Tree_Update_Topnode_Kicks();

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

#ifdef GRAVITY_TREE
		Float dv[3] = { dt*P[ipart].Acc[0], dt*P[ipart].Acc[1], 
			dt*P[ipart].Acc[2] };

		Gravity_Tree_Update_Kicks(dv, P[ipart].Tree_Parent); // kick tree nodes
#endif // GRAVITY_TREE
	}
	
	Gravity_Tree_Update_Topnode_Kicks();

	Profile("Second Kick");

	return ;
}
