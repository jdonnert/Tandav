#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "Gravity/gravity.h"

#ifdef GRAVITY
static void accel_gravity();
static inline void zero_active_particle_accelerations() {};
#else
static inline void accel_gravity() {};
static void zero_active_particle_accelerations();
#endif


/* 
 * Collect all accelerations on particles 
 */

void Compute_Acceleration()
{
	Profile("Accelerations");

	zero_active_particle_accelerations(); // ! GRAVITY

	accel_gravity(); // GRAVITY, needs previous P.Acc and overwrites it

	Profile("Accelerations");

	return ;
}

#ifdef GRAVITY

static void accel_gravity_tree()
{
	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	Gravity_Tree_Acceleration();

	if (Sig.First_Step)
		Gravity_Tree_Acceleration();

	Gravity_Tree_Periodic();

	return ;
}

static void accel_gravity_multi_grid()
{
	if (Sig.Fullstep)
		Gravity_Multi_Grid();

	return ;
}

static void accel_gravity()
{
	Profile("Gravity");

	accel_gravity_multi_grid(); // GRAVITY_MULTI_GRID

	accel_gravity_tree(); // GRAVITY_TREE

	Gravity_Simple_Accel(); // GRAVITY_SIMPLE, performs force test

	Profile("Gravity");

	return ;
}

#else // !GRAVITY

static void zero_active_particle_accelerations()
{
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;
	}

	return ;
}

#endif // ! GRAVITY
