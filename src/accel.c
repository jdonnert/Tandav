#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"

#ifdef GRAVITY
static void accel_gravity();
#endif

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
#endif

	Profile("Accelerations");

	return ;
}

#ifdef GRAVITY

static void accel_gravity()
{
	Profile("Gravity");

	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	Gravity_Tree_Acceleration();

	if (Sig.First_Step)
		Gravity_Tree_Acceleration();

	Gravity_Tree_Periodic();

	Gravity_Simple_Accel();

	if (Sig.Fullstep)
		Gravity_Multi_Grid_Long_Range();

	Profile("Gravity");

	return ;
}
#endif // ! GRAVITY

static void zero_active_particle_accelerations()
{
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;
	}

	return ;
}

