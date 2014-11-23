#include "../globals.h"
#include "gravity.h"

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

void Accel_Gravity_Simple()
{
	Profile("Gravity_Simple");

	rprintf("Direct Gravity ... ");

	#pragma omp for
	for (int i = 0; i < Task.Npart_Total; i++) {

		int ipart = Active_Particle_List[i];

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

	} // for ipart

	rprintf("done\n");
	
	Profile("Gravity_Simple");

	return ;
}
