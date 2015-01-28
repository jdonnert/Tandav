#include "../globals.h"
#include "gravity.h"

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static double mean_error = 0, max_error = 0; 
static int worst_part = -1;

/*
 * This computes the gravitational interaction via direct summation and shows
 * the relative error resp. the old force. Note that the max relative error can 
 * become large if one component of the force is close to 0 without consequence.
 * Check to the total force to make sure this is not the case.
 */

void Accel_Gravity_Simple()
{
	Profile("Gravity_Simple");

	rprintf("Direct Gravity ... ");

	mean_error = max_error = 0;
	worst_part = -1;

	#pragma omp for reduction(+:mean_error)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		double acc[3] = { P[ipart].Acc[0], P[ipart].Acc[1], P[ipart].Acc[2] };

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
		P[ipart].Grav_Acc[0] = P[ipart].Grav_Acc[1] = P[ipart].Grav_Acc[2] = 0;
#endif

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = 0;
#endif

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
			
				rinv = sqrt(14*u-84*u3+140*u2*u2-90*u2*u3+21*u3*u3)/h;
			} 
	
			double acc_mag = Const.Gravity * P[jpart].Mass * p2(rinv);

			P[ipart].Acc[0] += -acc_mag * dx * rinv;
			P[ipart].Acc[1] += -acc_mag * dy * rinv;
			P[ipart].Acc[2] += -acc_mag * dz * rinv;

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
			P[ipart].Grav_Acc[0] = -acc_mag * dx * rinv;
			P[ipart].Grav_Acc[1] = -acc_mag * dy * rinv;
			P[ipart].Grav_Acc[2] = -acc_mag * dz * rinv;
#endif

#ifdef GRAVITY_POTENTIAL
			if (r < h) {// WC2 kernel softening

				double u = r/h;
				double u2 = u*u;
				double u3 = u2*u;

				rinv = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/h;
			}

			P[ipart].Grav_Pot += -Const.Gravity * P[jpart].Mass *rinv;
#endif
		} // for jpart

		double error[3] = {(acc[0] - P[ipart].Acc[0]) / P[ipart].Acc[0],
							(acc[1] - P[ipart].Acc[1]) / P[ipart].Acc[1],
							(acc[2] - P[ipart].Acc[2]) / P[ipart].Acc[2] };
		
		double errorl = ALENGTH3(error);

		mean_error += errorl;

		#pragma omp critical
		{
		
		if (errorl > max_error) {

			max_error = ALENGTH3(error);
			worst_part = ipart;
		}

		} // omp critical
	
		//printf("ipart = %d %g | %g %g %g | %g %g %g \n", 
		//		ipart, errorl, acc[0], acc[1], acc[2],
		//		P[ipart].Acc[0],P[ipart].Acc[1],P[ipart].Acc[2] );

	} // for ipart

	rprintf("done\n");
	
	rprintf("\nForce test: max error %g @ %d, mean error %g \n\n", 
			max_error, worst_part, mean_error/NActive_Particles);

	Profile("Gravity_Simple");

	return ;
}
