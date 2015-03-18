#include "../globals.h"

#ifdef GRAVITY_SIMPLE 

#include "gravity.h"
#include "gravity_periodic.h"

static const double H = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static double Mean_Error = 0, Max_Error = 0;
static int Worst_Part = -1;

/*
 * This computes the gravitational interaction via direct summation and shows
 * the relative error resp. the old force. Note that the max relative error can 
 * become large if one component of the force is close to 0 without consequence.
 * Check to the total force to make sure this is not the case. Make sure the
 * mean rel. error in computations with many particles is a few percent if all
 * particles are kicked !
 */

void Gravity_Simple_Accel()
{
	Profile("Gravity_Simple");

	rprintf("Direct Gravity, get a coffee ... ");

	#pragma omp single
	{
	
	Mean_Error = Max_Error = 0;
	
	Worst_Part = -1;
	
	} // omp single

	#pragma omp for reduction(+:Mean_Error)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		double acc[3] = { P[ipart].Acc[0], P[ipart].Acc[1], P[ipart].Acc[2] };

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = 0;
#endif

		for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

			if (jpart == ipart)
				continue;

			Float dr[3] = { P[ipart].Pos[0] - P[jpart].Pos[0],
							P[ipart].Pos[1] - P[jpart].Pos[1],
							P[ipart].Pos[2] - P[jpart].Pos[2]};

			Float acc_periodic[3] = { 0 };

			Ewald_Correction(dr, &acc_periodic[0]); // PERIODIC

			P[ipart].Acc[0] += acc_periodic[0];
			P[ipart].Acc[1] += acc_periodic[1];
			P[ipart].Acc[2] += acc_periodic[2];

#ifdef GRAVITY_POTENTIAL
			Float pot_periodic = 0; 

			Ewald_Potential(dr, &pot_periodic); // PERIODIC

			P[ipart].Grav_Pot += pot_periodic;
#endif

			Periodic_Nearest(dr); // PERIODIC 

			double r = ALENGTH3(dr);

			double rinv = 1/r;

			if (r < H) {

				double u = r/H;
				double u2 = u*u;
				double u3 = u2*u;

				rinv = sqrt(14*u-84*u3+140*u2*u2-90*u2*u3+21*u3*u3)/H;
			} 

			double acc_mag = Const.Gravity * P[jpart].Mass * p2(rinv);

			P[ipart].Acc[0] += -acc_mag * dr[0] * rinv;
			P[ipart].Acc[1] += -acc_mag * dr[1] * rinv;
			P[ipart].Acc[2] += -acc_mag * dr[2] * rinv;
			
#ifdef GRAVITY_POTENTIAL
			if (r < H) { // WC2 kernel softening

				double u = r/H;
				double u2 = u*u;
				double u3 = u2*u;

				rinv = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/H;
			}

			P[ipart].Grav_Pot += -Const.Gravity * P[jpart].Mass *rinv;
#endif
		} // for jpart

		double error[3] = {(acc[0] - P[ipart].Acc[0]) / P[ipart].Acc[0],
							(acc[1] - P[ipart].Acc[1]) / P[ipart].Acc[1],
							(acc[2] - P[ipart].Acc[2]) / P[ipart].Acc[2] };
		
		double errorl = ALENGTH3(error);

		Mean_Error += errorl;

		#pragma omp critical
		{
		
		if (errorl > Max_Error) {

			Max_Error = ALENGTH3(error);
			Worst_Part = ipart;
		}

		} // omp critical
	
		//printf("ipart = %d %g | %g %g %g | %g %g %g \n", 
		//		ipart, errorl, acc[0], acc[1], acc[2],
		//		P[ipart].Acc[0],P[ipart].Acc[1],P[ipart].Acc[2] );

	} // for ipart

	rprintf("done\n");
	
	rprintf("\nForce test: NActive %d, max error %g @ %d, mean error %g \n\n", 
			NActive_Particles, Max_Error, Worst_Part, 
			Mean_Error/NActive_Particles);

	Profile("Gravity_Simple");

	return ;
}

#endif // GRAVITY_SIMPLE
