/* Collect all forces on particle ipart */
#include "globals.h"
#include "force.h"

static void Force_Gravity_Simple(const int ipart, float force[3]);

void Total_Force(const int ipart, float total_force[3])
{
	total_force[0] = total_force[1] = total_force[2] = 0;

	float force[3] = { 0 };

#ifdef GRAVITY
	Force_Gravity_Simple(ipart, force);
#endif // GRAVITY

	total_force[0] += force[0];
	total_force[1] += force[1];
	total_force[2] += force[2];

#ifdef OUTPUT_FORCE
	P[ipart].Force[0] = total_force[0];
	P[ipart].Force[1] = total_force[1];
	P[ipart].Force[2] = total_force[2];
#endif // OUTPUT_FORCE

	return ;
}

static void Force_Gravity_Simple(const int ipart, float *force)
{
	const double m_i = P[ipart].Mass;

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

	return ;
}
