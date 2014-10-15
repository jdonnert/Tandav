/* 
 * Collect all accelerations on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"

static void accel_gravity_simple(const int ipart, double force[3], 
		double *potential);

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

		accel_gravity_simple(ipart, grav_accel, &grav_potential);

		accel[0] += grav_accel[0];
		accel[1] += grav_accel[1];
		accel[2] += grav_accel[2];

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = grav_potential;
#endif // GRAVITY_POTENTIAL
#endif // GRAVITY
	
		P[ipart].Acc[0] = (Float) accel[0];
		P[ipart].Acc[1] = (Float) accel[1];
		P[ipart].Acc[2] = (Float) accel[2]; 
	}
	
	Profile("Accelerations");

	#pragma omp barrier

	return ;
}

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static void accel_gravity_simple(const int ipart, double *force, 
		double *potential)
{
	for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

		if (jpart == ipart)
			continue;

		double dx = P[ipart].Pos[0] - P[jpart].Pos[0];
		double dy = P[ipart].Pos[1] - P[jpart].Pos[1];
		double dz = P[ipart].Pos[2] - P[jpart].Pos[2];

		double r = sqrt(dx*dx + dy*dy + dz*dz);

		double rinv = 1/r;

		if (r < h) {
	
			double u = r/h;
			double u2 = u*u;
			double u3 = u2*u;
			
			rinv = sqrt(14*u- 84*u3 + 140 * u2*u2 - 90*u2*u3 + 21*u3*u3 )/h;
		} 
	
		double fmag = Const.Gravity * P[jpart].Mass * p2(rinv);

		force[0] += -fmag * dx * rinv;
		force[1] += -fmag * dy * rinv;
		force[2] += -fmag * dz * rinv;

#ifdef GRAVITY_POTENTIAL
		if (r < h) {// WC2 kernel softening

			double u = r/h;
			double u2 = u*u;
			double u3 = u2*u;

			rinv = (7*u2 - 21*u2*u2 + 28*u3*u2 - 15*u3*u3 + u3*u3*u*8 - 3)/h;
		}

		*potential += -Const.Gravity * P[jpart].Mass *rinv;
#endif
	}

	return ;
}
