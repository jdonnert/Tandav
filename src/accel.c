/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"

static void accel_gravity();

static void zero_active_particle_accelerations();

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

#ifdef GRAVITY_SIMPLE
	Accel_Gravity_Simple();
#endif

#ifdef GRAVITY_TREE
	if (Sig.Domain_Updated)
		Gravity_Tree_Build();

	Gravity_Tree_Acceleration();

#ifdef PERIODIC
	Gravity_Tree_Periodic();
#endif	

#endif // GRAVITY_TREE

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
