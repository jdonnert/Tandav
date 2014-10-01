/* 
 * Collect all forces on particle ipart 
 */

#include "globals.h"
#include "accel.h"
#include "timestep.h"

static void accel_gravity_simple(const int ipart, Float force[3], 
		Float *potential);

void Compute_Acceleration()
{
	Profile("Accelerations");

	rprintf("Computing Accelerations... ");
	
	#pragma omp parallel for 
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];
	
		Float accel[3] = { 0 };
		Float potential = 0;

#ifdef GRAVITY
		accel_gravity_simple(ipart, accel, &potential);
#endif // GRAVITY
	
		P[ipart].Acc[0] = accel[0];
		P[ipart].Acc[1] = accel[1];
		P[ipart].Acc[2] = accel[2]; 

		P[ipart].Potential = potential;
	}
	
	rprintf("done\n");

	Profile("Accelerations");

	return ;
}

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static void accel_gravity_simple(const int ipart, Float *force, 
		Float *potential)
{
	for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

		if (jpart == ipart)
			continue;

		double dx = P[ipart].Pos[0] - P[jpart].Pos[0];
		double dy = P[ipart].Pos[1] - P[jpart].Pos[1];
		double dz = P[ipart].Pos[2] - P[jpart].Pos[2];

		//double r =  sqrt(dx*dx + dy*dy + dz*dz + p2(GRAV_SOFTENING) );
	
		double r = sqrt(dx*dx + dy*dy + dz*dz);

		double rinv = 1/r;

		if (r < h) {
	
			double u = r/h;
			
			rinv = sqrt(u*(14 + u*u * (-84 + u * (140 + u * (-90 + 21*u)))))/h;
		} 
	
		double fmag = Const.Gravity * P[jpart].Mass * p2(rinv);

		force[0] += -fmag * dx * rinv;
		force[1] += -fmag * dy * rinv;
		force[2] += -fmag * dz * rinv;

		//r = sqrt(dx*dx + dy*dy + dz*dz);

	//printf("%g %g %g %g \n", 
//			sqrt(r2), rinv, fmag, 
//			sqrt(p2(force[0]) + p2(force[1]) + p2(force[2])));

	/*	if (r < h) { // WC2 kernel softening

			
			fmag *= (14*u3 - 84*u2*u3 + 140*u3*u3 - 90*u*u3*u3 + 21*u2*u3*u3);
			
			phi *= (7*u3 - 21*u2*u2 + 28*u3*u3 - 15*u3*u3*u + 8*u3*u3*u2 - 3*u);
		}
	*/	
		
		*potential += Const.Gravity * P[jpart].Mass *rinv;
	}

	return ;
}
