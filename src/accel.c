/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"

static void Accel_Gravity();

void Compute_Acceleration()
{
	Profile("Accelerations");

#ifdef GRAVITY
	Accel_Gravity();
#else
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) { 
			
		int ipart = Active_Particle_List[i];
	
		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;
	}
#endif // GRAVITY

	#pragma omp barrier

	Profile("Accelerations");

	return ;
}

static void Accel_Gravity()
{
	Profile("Gravity");

#ifdef GRAVITY_SIMPLE
	Accel_Gravity_Simple();
#endif

#ifdef GRAVITY_TREE
#pragma omp single
	Build_Gravity_Tree();
	
	Gravity_Tree_Acceleration();
#pragma omp barrier
#endif

	Profile("Gravity");

	return ;
}

