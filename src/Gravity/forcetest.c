#include "forcetest.h"

#ifdef GRAVITY_FORCETEST

#define NTEST 10

//static const Float H = 105/32 * GRAV_SOFTENING; // Plummer equiv softening

static Float H = 0;

static double Mean_Error = 0, Max_Error = 0;
static int Worst_Part = -1;

static double xacc_i = 0, yacc_i = 0, zacc_i = 0, pot_i = 0;
static int idx[NTEST] = { 0 };

/*
 * This computes the gravitational interaction via direct summation and shows
 * the relative error resp. the old force. Note that the total relative error 
 * can become large if one component of the force is close to 0 without 
 * consequence. We check NTEST particles only.
 */

void Gravity_Forcetest()
{
	Profile("Gravity_Forcetest");

	rprintf("Direct Gravity, N = %d \n"
			"      ID          Errors                 Tree Accel.         "
			"             Direct Accel. \n", NTEST);

	#pragma omp single
	{

	Mean_Error = Max_Error = 0;

	Worst_Part = -1;

	H = -41.0/32.0 * Param.Grav_Softening[1];

	} // omp single

	
	#pragma omp single
	for (int i = 0; i < NTEST; i++) 
		idx[i] = erand48(Task.Seed) * NActive_Particles;		

	for (int i = 0; i < NTEST; i++) {

		int ipart = Active_Particle_List[idx[i]];

		const double acc_old[3] = { P.Acc[0][ipart], 
									P.Acc[1][ipart], 
									P.Acc[2][ipart] };

		const double pos_i[3] = { P.Pos[0][ipart], 
								  P.Pos[1][ipart], 
								  P.Pos[2][ipart] };

		#pragma omp single 
		xacc_i = yacc_i = zacc_i = pot_i = 0;

		#pragma omp for reduction(+:xacc_i,yacc_i,zacc_i,pot_i)
		for (int jpart = 0; jpart < Sim.Npart_Total; jpart++) {

			Float dr[3] = { P.Pos[0][jpart] - pos_i[0],
							P.Pos[1][jpart] - pos_i[1],
							P.Pos[2][jpart] - pos_i[2] };

			Periodic_Nearest(dr); // PERIODIC 

			double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

			if (r2 > 0) {

				Float fac = Const.Gravity * P.Mass[jpart];
				Float fac_pot = Const.Gravity * P.Mass[jpart];

				if (r2 < Epsilon2[1]) { 

					Float u2 = r2 / Epsilon2[1];
	
					fac *= (175 - u2 * (294 - u2 * 135)) / (16*Epsilon3[1]) ;

					fac_pot *= (u2 * (175 - (u2 * 147  - u2 * 45)) - 105)
								/(32*Epsilon[1]);

				} else {

					Float r_inv = 1/SQRT(r2); // tempt the compiler for rsqrts

					fac *= r_inv * r_inv * r_inv;
					fac_pot *= r_inv;
				}

				xacc_i += fac * dr[0];
				yacc_i += fac * dr[1];
				zacc_i += fac * dr[2];

#ifdef PERIODIC
				Float result[3] = { 0 };
	
				Ewald_Correction(dr, &result[0]);

				xacc_i += Const.Gravity * P.Mass[jpart] * result[0];
				yacc_i += Const.Gravity * P.Mass[jpart] * result[1];
				zacc_i += Const.Gravity * P.Mass[jpart] * result[2];

#endif // PERIODIC

#ifdef GRAVITY_POTENTIAL
				if (r < H) { // WC2 kernel softening

					double u = r/H;
					double u2 = u*u;
					double u3 = u2*u;

					rinv = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/H;
				}

				pot_i += Const.Gravity * P.Mass[jpart] *rinv;
#endif

#if defined (GRAVITY_POTENTIAL) && defined(PERIODIC)
				Float pot_corr = 0;

				Ewald_Potential(dr, &pot_corr); // PERIODIC

				pot_i += Const.Gravity * P.Mass[jpart] * pot_corr;
#endif 

			} // r2 > 0
		} // for jpart

		#pragma omp single
		{

		double error[3] = { (acc_old[0] - xacc_i) / xacc_i,
							(acc_old[1] - yacc_i) / yacc_i,
							(acc_old[2] - zacc_i) / zacc_i };

		double error1 = sqrt( p2(error[0]) + p2(error[1]) + p2(error[2]));
	
		double acc = sqrt( p2(xacc_i) + p2(yacc_i) + p2(zacc_i) );
		double acc2 = sqrt( p2(acc_old[0]) + p2(acc_old[1]) + p2(acc_old[2]) );

		double error2 = (acc2 - acc) /acc;

		Mean_Error += error2;

		if (error2 > Max_Error) {

			Max_Error = error2;
			Worst_Part = ipart;
		}

		printf("   %9d  %05.4f %+05.4f   %+09.3f %+09.3f %+09.3f	"
				"%+09.3f %+09.3f %+09.3f\n",
				P.ID[ipart], error1, error2, acc_old[0], acc_old[1], 
				acc_old[2],	xacc_i, yacc_i, zacc_i);

		} // omp single
	} // for ipart

	rprintf("done\n");

	rprintf("\nForce test: NActive %d, max error %g @ %d, mean error %g \n\n",
			NActive_Particles, Max_Error, Worst_Part,
			Mean_Error/NTEST);

	Profile("Gravity_Forcetest");

	return ;
}

#endif // GRAVITY_FORCETEST


// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
