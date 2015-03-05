#include "globals.h"

#ifdef PERIODIC

void Periodic_Constrain_Particles_To_Box()
{
	const Float boxsize[3] = {Sim.Boxsize[0], Sim.Boxsize[1], Sim.Boxsize[2]};

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		while (P[ipart].Pos[0] < 0)
			P[ipart].Pos[0] += boxsize[0];

		while (P[ipart].Pos[0] >= boxsize[0])
			P[ipart].Pos[0] -= boxsize[0];

		while (P[ipart].Pos[1] < 0)
			P[ipart].Pos[1] += boxsize[1];

		while (P[ipart].Pos[1] >= boxsize[1])
			P[ipart].Pos[1] -= boxsize[1];

		while (P[ipart].Pos[2] < 0)
			P[ipart].Pos[2] += boxsize[2];

		while (P[ipart].Pos[2] >= boxsize[2])
			P[ipart].Pos[2] -= boxsize[2];
	} // for i

	return ;
}

/*
 * Map a distance "dx" on component "i" to the nearest distance given 
 * the periodic box. This is not very fast.
 */

Float Periodic_Nearest(const Float dx)
{
	Float dx_periodic = dx;

	if (dx > 0.5 * Sim.Boxsize[0])
		dx_periodic = dx - Sim.Boxsize[0];
	else if (dx < -0.5 * Sim.Boxsize[0])
		dx_periodic = dx + Sim.Boxsize[0];

	return dx_periodic;
}

Float Periodic_Nearest_Noncubic(const Float dx, const int i)
{
	Float dx_periodic = dx;

	if (dx > 0.5 * Sim.Boxsize[i])
		dx_periodic = dx - Sim.Boxsize[i];
	else if (dx < -0.5 * Sim.Boxsize[i])
		dx_periodic = dx + Sim.Boxsize[i];

	return dx_periodic;
}


#endif // PERIODIC
