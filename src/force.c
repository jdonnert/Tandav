/* Collect all forces on particle ipart */
#include "globals.h"
#include "force.h"
#include "timestep.h"

static void force_Gravity_Simple(const int ipart, Float force[3], 
		Float *potential);

void Compute_Forces()
{
	Profile("Forces");

	for (int ipart = 0; ipart < Task.NpartTotal; ipart++){
	
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

const double h = pow( 21/(2*PI*GRAV_SOFTENING) ,1/3); // Plummer

static void force_Gravity_Simple(const int ipart, Float *force, 
		Float *potential)
{
	const double m_i = P[ipart].Mass;

	Profile("Gravity Simple");
	
	for (int jpart = 0; jpart < Task.NpartTotal; jpart++ ) {
	
		if (jpart == ipart)
			continue;

		double mpart = P[jpart].Mass;
		
		double dx = P[ipart].Pos[0] - P[jpart].Pos[0];
		double dy = P[ipart].Pos[1] - P[jpart].Pos[1];
		double dz = P[ipart].Pos[2] - P[jpart].Pos[2];

		double r = sqrt( dx*dx + dy*dy + dz*dz );

		double fmag = Const.Gravity * mpart * m_i/ p2(r);

		if (r < h) {
		
			fmag = 1;


		}


		force[0] += -fmag * dx/r;
		force[1] += -fmag * dy/r;
		force[2] += -fmag * dz/r;
	}
	
	Profile("Gravity Simple");

	return ;
}

static double softening_kernel_force(const double r, const double h)
{
	return 1;
}

static double softening_kernel_potential(const double r, const double h)
{
	return 1;
}
