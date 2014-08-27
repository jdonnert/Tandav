/* Collect all forces on particle ipart */
#include "globals.h"
#include "force.h"
#include "timestep.h"

static void force_Gravity_Simple(const int ipart, Float force[3], 
		Float *potential);

void Compute_Forces()
{
	Profile("Forces");

	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {
	
		Float force[3] = { 0 };
		Float potential = 0;

#ifdef GRAVITY
		force_Gravity_Simple(ipart, force, &potential);
#endif // GRAVITY
	
		P[ipart].Force[0] = force[0];
		P[ipart].Force[1] = force[1];
		P[ipart].Force[2] = force[2];

		P[ipart].Potential = potential;
	}

	Profile("Forces");

	return ;
}

const double h = 3 / GRAV_SOFTENING; // Plummer equivalent softening

static void force_Gravity_Simple(const int ipart, Float *force, 
		Float *potential)
{
	const double m_i = P[ipart].Mass;
	
	Profile("Gravity Simple");
	
	for (int jpart = 0; jpart < Task.NpartTotal; jpart++) {
	
		if (jpart == ipart)
			continue;

		double mpart = P[jpart].Mass;
		
		double dx = P[ipart].Pos[0] - P[jpart].Pos[0];
		double dy = P[ipart].Pos[1] - P[jpart].Pos[1];
		double dz = P[ipart].Pos[2] - P[jpart].Pos[2];

		double r = sqrt( dx*dx + dy*dy + dz*dz );

		double fmag = Const.Gravity * mpart * m_i / p2(r);
		double phi = Const.Gravity * mpart / r;

		if (r < h) { // WC2 kernel softening

			double u = r/h;
			double u2 = u*u;
			double u3 = u2*u;

			fmag *= (14*u3 - 84*u2*u3 + 140*u3*u3 - 90*u*u3*u3 + 21*u2*u3*u3);
			
			phi *= (7*u3 - 21*u2*u2 + 28*u3*u3 - 15*u3*u3*u + 8*u3*u3*u2 - 3*u);
		}

		force[0] += -fmag * dx/r;
		force[1] += -fmag * dy/r;
		force[2] += -fmag * dz/r;

		*potential = phi;
	}
	
	Profile("Gravity Simple");

	return ;
}
