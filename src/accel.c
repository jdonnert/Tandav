#include "accel.h"

#ifdef GRAVITY
static void accel_gravity();
static inline void zero_active_particle_accelerations() {};
#else
static inline void gravity_accel() {};
static void zero_active_particle_accelerations();
#endif

/* 
 * Collect all accelerations on particles. Note that GRAVITY needs the previous
 * accelerations and hence we do not zero P.Accel if it is switched on.
 */

void Compute_Acceleration()
{
	zero_active_particle_accelerations(); // !GRAVITY

	accel_gravity(); // GRAVITY

	// accel_hydro();

	return ;
}

#ifdef GRAVITY

static void gravity_accel_tree()
{
	if (Sig.Tree_Update)
		Gravity_Tree_Build();

	if (Sig.Prepare_Step) {

		Sig.Use_BH_Criterion = true;
		
		Gravity_Tree_Acceleration();
	
		Sig.Use_BH_Criterion = false;
	}

	Gravity_Tree_Acceleration();

	return ;
}

static void accel_gravity()
{
	gravity_accel_tree(); // GRAVITY_TREE
	
	//gravity_accel_fmm(); // GRAVITY_FMM
	
	//gravity_accel_mg();

	Gravity_Simple_Accel(); // GRAVITY_SIMPLE, performs force test

	return ;
}

#endif 

#ifndef GRAVITY
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
