#include "accel.h"

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
	Safe_Last_Accel(); // and zero P.Acc

	Gravity_Tree_Acceleration(); // GRAVITY_TREE

	Gravity_FMM_Acceleration(); // GRAVITY_FMM

	Gravity_Forcetest(); // GRAVITY_FORCETEST

	// Hydro_Acceleration(); // HYDRO

	return ;
}

void Safe_Last_Accel()
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






Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
