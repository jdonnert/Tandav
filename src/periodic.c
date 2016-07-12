#include "periodic.h"

#ifdef PERIODIC

static Float Boxhalf[3] = { 0 };

void Periodic_Setup()
{
#ifndef PERIODIC_NO_CUBE
	
	Boxhalf[0] = Boxhalf[1] = Boxhalf[2] = Sim.Boxsize[0]/2.0;

	rprintf("Periodic Boxsize = %g, Boxhalf = %g \n\n", 
			Sim.Boxsize[0], Boxhalf[0]);

#else // PERIODIC_NO_CUBE
	
	Boxhalf[0] = Sim.Boxsize[0];
	Boxhalf[1] = Sim.Boxsize[1];
	Boxhalf[2] = Sim.Boxsize[2];
	
	rprintf("Periodic Boxsizes = %g %g %g \n\n", Sim.Boxsize[0], 
			Sim.Boxsize[1],Sim.Boxsize[2]);

#endif // PERIODIC_NO_CUBE

	#pragma omp parallel
	Periodic_Constrain_Particles_To_Box(); // PERIODIC

	return ;
}

void Periodic_Constrain_Particles_To_Box()
{
	#pragma omp for 
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		while (P.Pos[0][ipart] < 0)
			P.Pos[0][ipart] += Sim.Boxsize[0];

		while (P.Pos[0][ipart] >= Sim.Boxsize[0])
			P.Pos[0][ipart] -= Sim.Boxsize[0];

#ifdef PERIODIC_NO_CUBE
		while (P.Pos[1][ipart] < 0)
			P.Pos[1][ipart] += Sim.Boxsize[1];

		while (P.Pos[1][ipart] >= Sim.Boxsize[1])
			P.Pos[1][ipart] -= Sim.Boxsize[1];

		while (P.Pos[2][ipart] < 0)
			P.Pos[2][ipart] += Sim.Boxsize[2];

		while (P.Pos[2][ipart] >= Sim.Boxsize[2])
			P.Pos[2][ipart] -= Sim.Boxsize[2];
#endif
	} // for i

	return ;
}

/*
 * Do the periodic mapping on a 3D distance array. We are relying on link time 
 * optimization of the compiler to do the inlining & optimization for us. 
 */

void Periodic_Nearest(Float dr[3])
{
	if (dr[0] > Boxhalf[0])
		dr[0] -= Sim.Boxsize[0];
	else if (dr[0] < -Boxhalf[0])
		dr[0] += Sim.Boxsize[0];

	if (dr[1] > Boxhalf[1])
		dr[1] -= Sim.Boxsize[1];
	else if (dr[1] < -Boxhalf[1])
		dr[1] += Sim.Boxsize[1];

	if (dr[2] > Boxhalf[2])
		dr[2] -= Sim.Boxsize[2];
	else if (dr[2] < -Boxhalf[2])
		dr[2] += Sim.Boxsize[2];

	return ;
}
#endif // PERIODIC

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
