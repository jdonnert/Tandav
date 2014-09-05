/* Collect all forces on particle ipart */
#include "globals.h"
#include "force.h"
#include "timestep.h"

static void force_Gravity_Simple(const int ipart, Float force[3], 
		Float *potential);

void Compute_Forces()
{
	Profile("Forces");

	rprintf("Computing Forces ... ");
	
#pragma omp parallel for 
	for (int i = 0; i < NActiveParticles; i++) {
		
		int ipart = ActiveParticleList[i];
		
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
	
	rprintf("done\n");

	Profile("Forces");

	return ;
}

const static double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static void force_Gravity_Simple(const int ipart, Float *force, 
		Float *potential)
{
	for (int jpart = 0; jpart < Sim.NpartTotal; jpart++) {
	
		if (jpart == ipart)
			continue;

		double mpart = P[jpart].Mass;
		
		double dx = P[ipart].Pos[0] - P[jpart].Pos[0];
		double dy = P[ipart].Pos[1] - P[jpart].Pos[1];
		double dz = P[ipart].Pos[2] - P[jpart].Pos[2];

		double r2 = dx*dx + dy*dy + dz*dz;

		double r2inv = 1 / r2;

		if (r2 < h*h) {
		
			double u = sqrt(r2)/h;
			double u2 = u*u;
			double u3 = u2*u;
			
			r2inv = (14*u - 84*u3 + 140*u2*u2 - 90*u3*u2 + 21*u3*u3 ) / h /h;
		} 

		double fmag = Const.Gravity * mpart * P[ipart].Mass * r2inv;

		double rinv = sqrt(r2inv);

		force[0] += -fmag * dx* rinv;
		force[1] += -fmag * dy* rinv;
		force[2] += -fmag * dz* rinv;

	//printf("%g %g %g %g \n", 
//			sqrt(r2), rinv, fmag, 
//			sqrt(p2(force[0]) + p2(force[1]) + p2(force[2])));

	/*	if (r < h) { // WC2 kernel softening

			

			fmag *= (14*u3 - 84*u2*u3 + 140*u3*u3 - 90*u*u3*u3 + 21*u2*u3*u3);
			
			phi *= (7*u3 - 21*u2*u2 + 28*u3*u3 - 15*u3*u3*u + 8*u3*u3*u2 - 3*u);
		}
	*/	
		double phi = Const.Gravity * mpart * rinv;
		
		*potential += phi;
	}

	return ;
}
