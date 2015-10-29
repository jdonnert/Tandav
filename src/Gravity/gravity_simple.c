#include "../globals.h"

#ifdef GRAVITY_SIMPLE 

#include "gravity.h"
#include "gravity_periodic.h"

//static const double H = GRAV_SOFTENING / 3.0; // Plummer equivalent softening
static const Float H = 105/32 * GRAV_SOFTENING; // Plummer equiv softening

static double Mean_Error = 0, Max_Error = 0;
static int Worst_Part = -1;

/*
 * This computes the gravitational interaction via direct summation and shows
 * the relative error resp. the old force. Note that the max relative error 
 * can become large if one component of the force is close to 0 without 
 * consequence.
 * Check to the total force to make sure this is not the case. Make sure the
 * mean rel. error in computations with many particles is a few percent if all
 * particles are kicked !
 */

void Gravity_Simple_Accel()
{
	Profile("Gravity_Simple");

	rprintf("Direct Gravity, get a coffee ... \n");

	#pragma omp single
	{

	Mean_Error = Max_Error = 0;

	Worst_Part = -1;

	} // omp single

	int cnt = 0;

	#pragma omp for reduction(+:Mean_Error)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		if (P[ipart].ID > 5)
			continue;

		cnt++;

		double acc[3] = { P[ipart].Acc[0], P[ipart].Acc[1], P[ipart].Acc[2] };

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = 0;
#endif

		for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

			double dr[3] = { P[jpart].Pos[0] - P[ipart].Pos[0],
							P[jpart].Pos[1] - P[ipart].Pos[1],
							P[jpart].Pos[2] - P[ipart].Pos[2]};

			Periodic_Nearest(dr); // PERIODIC 

			double r = ALENGTH3(dr);

			double rinv = 1/r;

			if (r < H) {

				double u = r/H;
				double u2 = u*u;

				rinv = sqrt(u * (135*u2*u2 - 294*u2 + 175))/(4*H) ;
				
				//double u3 = u2*u;
				//rinv = sqrt(14*u-84*u3+140*u2*u2-90*u2*u3+21*u3*u3)/H;
			}

			double acc_mag = Const.Gravity * P[jpart].Mass * p2(rinv);

			P[ipart].Acc[0] += acc_mag * dr[0] * rinv;
			P[ipart].Acc[1] += acc_mag * dr[1] * rinv;
			P[ipart].Acc[2] += acc_mag * dr[2] * rinv;

#ifdef PERIODIC
			Float ew_corr[3] = { 0 };

			Ewald_Correction(dr, &ew_corr[0]);

			P[ipart].Acc[0] += Const.Gravity * P[jpart].Mass * ew_corr[0];
			P[ipart].Acc[1] += Const.Gravity * P[jpart].Mass * ew_corr[1];
			P[ipart].Acc[2] += Const.Gravity * P[jpart].Mass * ew_corr[2];

#endif // PERIODIC

#ifdef GRAVITY_POTENTIAL
			if (r < H) { // WC2 kernel softening

				double u = r/H;
				double u2 = u*u;
				double u3 = u2*u;

				rinv = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/H;
			}

			P[ipart].Grav_Pot += Const.Gravity * P[jpart].Mass *rinv;
#endif

#if  defined (GRAVITY_POTENTIAL) && defined(PERIODIC)
			Float pot_corr = 0;

			Ewald_Potential(dr, &pot_corr); // PERIODIC

			P[ipart].Grav_Pot += Const.Gravity * P[jpart].Mass * pot_corr;
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

		printf("  ipart=%d ID=%d err=%g tree acc=%g %g %g dir acc=%g %g %g \n",
				ipart, P[ipart].ID, errorl, acc[0], acc[1], acc[2],
				P[ipart].Acc[0],P[ipart].Acc[1],P[ipart].Acc[2] );

	} // for ipart

	rprintf("done\n");

	rprintf("\nForce test: NActive %d, max error %g @ %d, mean error %g \n\n",
			NActive_Particles, Max_Error, Worst_Part,
			Mean_Error/cnt);

	Profile("Gravity_Simple");
	
	return ;
}

#endif // GRAVITY_SIMPLE
