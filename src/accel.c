#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"

static void accel_gravity();

static void zero_active_particle_accelerations();

/* 
 * Collect all accelerations on particles 
 */

void Compute_Acceleration()
{
	Profile("Accelerations");

#ifdef GRAVITY
	accel_gravity(); // needs previous P.Acc and overwrites it
#else
	zero_active_particle_accelerations();
#endif // GRAVITY

	Profile("Accelerations");

	return ;
}

static void accel_gravity()
{
	Profile("Gravity");

#ifdef GRAVITY_TREE

	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	Gravity_Tree_Acceleration();
	
	if (Sig.First_Step) 
		Gravity_Tree_Acceleration();

#ifdef PERIODIC
	 Gravity_Tree_Ewald_Correction();
#endif	

#endif // GRAVITY_TREE

#ifdef GRAVITY_SIMPLE
	Accel_Gravity_Simple();
#endif

#ifdef GRAVITY_GRID
	if (Sig.Full_Step)
		Gravity_Grid_Long_Range
#endif
	Profile("Gravity");

	return ;
}


static void zero_active_particle_accelerations()
{
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) { 
			
		int ipart = Active_Particle_List[i];
	
		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;
	}
		
	return ;
}
