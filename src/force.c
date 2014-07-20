/* Collect all forces on particle ipart */
#include "globals.h"
#include "force.h"
#include "timestep.h"

static void force_Gravity_Simple(const int ipart, float force[3]);

void Compute_Forces()
{
	Profile("Forces");

	for (int ipart = 0; ipart < Task.NpartTotal; ipart++){
	
		float force[3] = { 0 };
	
#ifdef GRAVITY
		force_Gravity_Simple(ipart, force);
#endif // GRAVITY
	
		float last_acc = len3(P[ipart].Force) / P[ipart].Mass;
		float acc =  len3(force) / P[ipart].Mass;

		P[ipart].Surge = (acc - last_acc) / Time.Step;
		
		P[ipart].Force[0] = force[0];
		P[ipart].Force[1] = force[1];
		P[ipart].Force[2] = force[2];
	}

	Profile("Forces");

	return ;
}

static void force_Gravity_Simple(const int ipart, float *force)
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

		force[0] += -fmag * dx/r;
		force[1] += -fmag * dy/r;
		force[2] += -fmag * dz/r;
	}
	
	Profile("Gravity Simple");

	return ;
}
