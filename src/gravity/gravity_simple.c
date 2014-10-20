#include "../globals.h"
#include "gravity.h"

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

void Accel_Gravity_Simple(const int ipart, double *accel, double *potential)
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
	
		double acc_mag = Const.Gravity * P[jpart].Mass * p2(rinv);

		accel[0] += -acc_mag * dx * rinv;
		accel[1] += -acc_mag * dy * rinv;
		accel[2] += -acc_mag * dz * rinv;

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
