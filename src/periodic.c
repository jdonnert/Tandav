#include "globals.h"

#ifdef PERIODIC

void Periodic_Constrain_Particles_To_Box()
{
	const Float boxsize[3] = {Sim.Boxsize[0], Sim.Boxsize[1], Sim.Boxsize[2]};

	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

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
 * Do the periodic mapping on a 3D distance array. We are relying on link time 
 * optimization of the compiler to do the inlining for us. Make sure to put
 * the appropriate compiler switches.
 */

void Periodic_Nearest(Float dr[3])
{
	for (int i = 0; i < 3; i++) {
	
		if (dr[i] > 0.5 * Sim.Boxsize[i])
			dr[i] -= Sim.Boxsize[i];
		else if (dr[i] < -0.5 * Sim.Boxsize[i])
			dr[i] += Sim.Boxsize[i];
	}

	return ;
}
#endif // PERIODIC
