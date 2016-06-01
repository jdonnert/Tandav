#include "accel.h"

static void zero_active_particle_accelerations();

#ifdef GRAVITY
static void accel_gravity();
#else
static inline void gravity_accel() {};
#endif // GRAVITY

/* 
 * Collect all accelerations on particles.
 */

void Compute_Acceleration()
{
	zero_active_particle_accelerations(); // also sets P.Last_Acc_Mag

	accel_gravity(); // GRAVITY

	Gravity_Forcetest(); // GRAVITY_FORCETEST

	// accel_hydro();

	return ;
}

static void zero_active_particle_accelerations()
{
	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
	
		P.Last_Acc_Mag[ipart] = SQRT( p2(P.Acc[0][ipart]) 
									+ p2(P.Acc[1][ipart]) 
									+ p2(P.Acc[2][ipart]));

		P.Acc[0][ipart] = P.Acc[1][ipart] = P.Acc[2][ipart] = 0;
	}

	return ;
}

/*
 * Gravity driving routines. Depending on the algorithm, we execute different
 * gravity modules or combinations of them.
 */

#if defined(GRAVITY) && defined(GRAVITY_TREE) // pure tree
static void accel_gravity() 
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
#endif  // GRAVITY && GRAVITY_TREE

#if defined(GRAVITY) && defined(GRAVITY_FMM) // pure FMM
static void accel_gravity() 
{
	//Gravity_FMM_Build(); 
	
	//Gravity_FMM_Accel();

	return ;
}
#endif  // GRAVITY && GRAVITY_FMM


