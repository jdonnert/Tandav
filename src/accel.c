/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"
#include "memory.h"
static void Accel_Gravity();

void Compute_Acceleration()
{
	Profile("Accelerations");

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
			
		int ipart = Active_Particle_List[i];
	
		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[1] = 0;
	}

#ifdef GRAVITY
	Accel_Gravity();
#endif // GRAVITY

	Profile("Accelerations");

	return ;
}

static void Accel_Gravity()
{

#ifdef GRAVITY_SIMPLE
	Accel_Gravity_Simple();
#endif

#ifdef GRAVITY_TREE
	Build_Gravity_Tree();

	Gravity_Tree_Acceleration();

#endif


	return ;
}

