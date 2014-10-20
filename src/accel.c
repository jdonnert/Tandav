/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"
#include "gravity/gravity.h"


void Compute_Acceleration()
{
	Profile("Accelerations");

	#pragma omp for 
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];
	
		double accel[3] = { 0 };

#ifdef GRAVITY
		double grav_accel[3] = { 0 };
		double grav_potential = 0;

#ifdef GRAVITY_SIMPLE
		Accel_Gravity_Simple(ipart, grav_accel, &grav_potential);
#endif

#ifdef GRAVITY_TREE
		Accel_Gravity_Tree(ipart, grav_accel, &grav_potential);
#endif

#ifdef GRAVITY_PM
		Accel_Gravity_PM(ipart, grav_accel, &grav_potential);
#endif

		accel[0] += grav_accel[0];
		accel[1] += grav_accel[1];
		accel[2] += grav_accel[2];

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = grav_potential;
#endif 

#endif // GRAVITY
	
		P[ipart].Acc[0] = (Float) accel[0];
		P[ipart].Acc[1] = (Float) accel[1];
		P[ipart].Acc[2] = (Float) accel[2]; 
	}
	
	Profile("Accelerations");

	return ;
}



