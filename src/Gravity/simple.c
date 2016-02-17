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

	#pragma omp for simd reduction(+:Mean_Error)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		if (P.ID[ipart] > 10)
			continue;

		cnt++;

		const double acc[3] = { P.Acc[0][ipart], P.Acc[1][ipart], 
			P.Acc[2][ipart] };
		const double pos_i[3] = { P.Pos[0][ipart], P.Pos[1][ipart], 
			P.Pos[2][ipart] };

		double acc_i[3] = { 0 }; 

#ifdef GRAVITY_POTENTIAL
		P.Grav_Pot[ipart] = 0;
#endif

		for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

			Float dr[3] = {P.Pos[0][jpart] - pos_i[0],
							P.Pos[1][jpart] - pos_i[1],
							P.Pos[2][jpart] - pos_i[2]};

			Periodic_Nearest(dr); // PERIODIC 

			double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

			if (r2 > 0) {

				double r = sqrt(r2);
				double rinv = 1/r;
	
				if (r < H) {

					double u = r/H;
					double u2 = u*u;

					rinv = sqrt(u * (135*u2*u2 - 294*u2 + 175))/(4*H) ;
				}

				double acc_mag = Const.Gravity * P.Mass[jpart] * p2(rinv);
			
				acc_i[0] += acc_mag * dr[0] * rinv;
				acc_i[1] += acc_mag * dr[1] * rinv;
				acc_i[2] += acc_mag * dr[2] * rinv;

#ifdef PERIODIC
				Float result[3] = { 0 };

				Ewald_Correction(dr, &result[0]);

				acc_i[0] += Const.Gravity * P.Mass[jpart] * result[0];
				acc_i[1] += Const.Gravity * P.Mass[jpart] * result[1];
				acc_i[2] += Const.Gravity * P.Mass[jpart] * result[2];

#endif // PERIODIC

#ifdef GRAVITY_POTENTIAL
				if (r < H) { // WC2 kernel softening

					double u = r/H;
					double u2 = u*u;
					double u3 = u2*u;

					rinv = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/H;
				}

				P.Grav_Pot[ipart] += Const.Gravity * P.Mass[jpart] *rinv;
#endif

#if defined (GRAVITY_POTENTIAL) && defined(PERIODIC)
				Float pot_corr = 0;

				Ewald_Potential(dr, &pot_corr); // PERIODIC

				P.Grav_Pot[ipart] += Const.Gravity * P.Mass[jpart] * pot_corr;
#endif 

			} // r2 > 0
		} // for jpart

		P.Acc[0][ipart] = acc_i[0];
		P.Acc[1][ipart] = acc_i[1];
		P.Acc[2][ipart] = acc_i[2];

		double error[3] = { (acc[0] - acc_i[0]) / acc_i[0],
							(acc[1] - acc_i[1]) / acc_i[1],
							(acc[2] - acc_i[2]) / acc_i[2] };

		double errorl = ALENGTH3(error);

		Mean_Error += errorl;

		#pragma omp critical
		{

		if (errorl > Max_Error) {

			Max_Error = ALENGTH3(error);
			Worst_Part = ipart;
		}

		} // omp critical

		printf("  ipart=%d ID=%d pos=%g %g %g err=%g tree acc=%g %g %g "
				"dir acc=%g %g %g \n",
				ipart, P.ID[ipart], P.Pos[0][ipart], P.Pos[1][ipart],
				P.Pos[2][ipart], errorl, acc[0], acc[1], acc[2],
				P.Acc[0][ipart], P.Acc[1][ipart], P.Acc[2][ipart] );

	} // for ipart

	rprintf("done\n");

	rprintf("\nForce test: NActive %d, max error %g @ %d, mean error %g \n\n",
			NActive_Particles, Max_Error, Worst_Part,
			Mean_Error/cnt);

	Profile("Gravity_Simple");

	exit(0);
	return ;
}

#endif // GRAVITY_SIMPLE
