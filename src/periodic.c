#include "globals.h"

#ifdef PERIODIC

Float Boxsize = 0, Boxhalf = 0;

void Init_Periodic()
{
	Boxsize = Sim.Boxsize[2] = Sim.Boxsize[1] = Sim.Boxsize[0];
	Boxhalf = Boxsize/2.0;
	
	rprintf("Periodic Boxsize = %g, Boxhalf = %g \n\n", Boxsize, Boxhalf);

	return ;
}

void Periodic_Constrain_Particles_To_Box()
{
	const Float boxsize[3] = {Sim.Boxsize[0], Sim.Boxsize[1], Sim.Boxsize[2]};

	#pragma omp for 
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		while (P.Pos[0][ipart] < 0)
			P.Pos[0][ipart] += boxsize[0];

		while (P.Pos[0][ipart] >= boxsize[0])
			P.Pos[0][ipart] -= boxsize[0];

		while (P.Pos[1][ipart] < 0)
			P.Pos[1][ipart] += boxsize[1];

		while (P.Pos[1][ipart] >= boxsize[1])
			P.Pos[1][ipart] -= boxsize[1];

		while (P.Pos[2][ipart] < 0)
			P.Pos[2][ipart] += boxsize[2];

		while (P.Pos[2][ipart] >= boxsize[2])
			P.Pos[2][ipart] -= boxsize[2];
	} // for i

	return ;
}

/*
 * Do the periodic mapping on a 3D distance array. We are relying on link time 
 * optimization of the compiler to do the inlining & optimization for us. 
 */

void Periodic_Nearest(Float dr[3])
{
	if (dr[0] > Boxhalf)
		dr[0] -= Boxsize;
	else if (dr[0] < -Boxhalf)
		dr[0] += Boxsize;

	if (dr[1] > Boxhalf)
		dr[1] -= Boxsize;
	else if (dr[1] < -Boxhalf)
		dr[1] += Boxsize;

	if (dr[2] > Boxhalf)
		dr[2] -= Boxsize;
	else if (dr[2] < -Boxhalf)
		dr[2] += Boxsize;

	return ;
}
#endif // PERIODIC
