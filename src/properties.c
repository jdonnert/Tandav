#include "globals.h"

static void find_center_of_mass(double CoM_out[3]);
static void find_total_mass(double mass_out[1]);
static void find_total_kinetic_energy(double Ekin_out[1]);
static void find_angular_momentum(double ang_p_out[3]);
static void find_momentum(double mom_out[3]);

void Compute_Global_Simulation_Properties()
{
	Profile("Properties");

	find_total_mass(&Sim.Total_Mass);

	find_center_of_mass(&Sim.Center_Of_Mass[0]);

	find_total_kinetic_energy(&Sim.Kinetic_Energy);

	find_angular_momentum(&Sim.Angular_Momentum[0]);
	
	find_momentum(&Sim.Momentum[0]);

	Profile("Properties");

	return ;
}

static double Total_Mass = 0;

static void find_total_mass(double mass_out[1])
{
	#pragma omp single
	Total_Mass = 0;

	#pragma omp for reduction(+:Total_Mass) 
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		Total_Mass += P.Mass[ipart];

	#pragma omp single
	{

	MPI_Allreduce(MPI_IN_PLACE, &Total_Mass, 1, MPI_DOUBLE, MPI_SUM,
				  MPI_COMM_WORLD);

	mass_out[0] = Total_Mass;

	} // omp single

	return ;
}

static double CoM_X = 0, CoM_Y = 0, CoM_Z = 0; // can't reduce on array

static void find_center_of_mass(double CoM_out[3])
{
	#pragma omp single
	CoM_X = CoM_Y = CoM_Z = 0;

	#pragma omp for reduction(+:CoM_X,CoM_Y,CoM_Z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		CoM_X += P.Mass[ipart] * P.Pos[0][ipart];
		CoM_Y += P.Mass[ipart] * P.Pos[1][ipart];
		CoM_Z += P.Mass[ipart] * P.Pos[2][ipart];
	}

	#pragma omp single
	{

	double global_com[3] = { CoM_X, CoM_Y, CoM_Z  };

	MPI_Allreduce(MPI_IN_PLACE, global_com, 3, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	CoM_out[0] = global_com[0] / Sim.Total_Mass;
	CoM_out[1] = global_com[1] / Sim.Total_Mass;
	CoM_out[2] = global_com[2] / Sim.Total_Mass;

	} // omp single

	return ;

}

static double Ekin = 0;

static void find_total_kinetic_energy(double Ekin_out[1])
{
	#pragma omp single
	Ekin = 0;

	#pragma omp for reduction(+:Ekin)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		double v2 = p2(P.Vel[0][ipart])+p2(P.Vel[1][ipart])+p2(P.Vel[2][ipart]);
		
		Ekin += 0.5 * P.Mass[ipart] * v2;
	}
	
	#pragma omp single
	{

	MPI_Allreduce(Ekin_out, &Ekin, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	} // omp single

	return ;
}

static double Ang_p_x = 0, Ang_p_y = 0, Ang_p_z = 0;

static void find_angular_momentum(double ang_p_out[3])
{
	#pragma omp single
	Ang_p_x = Ang_p_y = Ang_p_z = 0;

	#pragma omp for reduction(+:Ang_p_x, Ang_p_y, Ang_p_z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		Ang_p_x += P.Mass[ipart] * (P.Pos[1][ipart]*P.Vel[2][ipart]
				- P.Pos[2][ipart]*P.Vel[1][ipart]);
		Ang_p_y += P.Mass[ipart] * (P.Pos[2][ipart]*P.Vel[0][ipart]
				- P.Pos[0][ipart]*P.Vel[2][ipart]);
		Ang_p_z += P.Mass[ipart] * (P.Pos[0][ipart]*P.Vel[1][ipart]
				- P.Pos[1][ipart]*P.Vel[0][ipart]);
	}

	#pragma omp single
	{

	double global_ang_p[3] = { Ang_p_x, Ang_p_y, Ang_p_z };

	MPI_Allreduce(ang_p_out, global_ang_p, 3, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	} // omp single

	return ;
}
	
static double mom_x = 0, mom_y = 0, mom_z = 0;

static void find_momentum(double mom_out[3]) 
{
	#pragma omp single
	mom_x = mom_y = mom_z = 0;
	
	#pragma omp for reduction(+:mom_x, mom_y, mom_z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		mom_x += P.Mass[ipart] * P.Vel[0][ipart];
		mom_y += P.Mass[ipart] * P.Vel[0][ipart];
		mom_z += P.Mass[ipart] * P.Vel[0][ipart];
	}

	#pragma omp single
	{

	double global_mom[3] = { mom_x, mom_y, mom_z };

	MPI_Allreduce(mom_out, global_mom, 3, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	} // omp single

	return ;
}
