#include "globals.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "Gravity/gravity_periodic.h"

struct Parameters_From_File Param;
struct Global_Simulation_Properties Sim;

#pragma omp threadprivate(Task)
struct Local_Task_Properties Task = { 0 };

struct Particle_Data *P = NULL;

void Read_and_Init()
{
	Init_Profiler();

	Profile("Init");

	Read_Parameter_File(Param.File);

	Init_Memory_Management();

	Init_Logs();

	switch (Param.Start_Flag) {

	case 0:

		Read_Snapshot(Param.Input_File); // also init particle structures

		#pragma omp parallel
		Sig.First_Step = Sig.Fullstep = true;

		break;

	case 1:

		Read_Restart_File();

		break;

	default:

		Assert(0, "Start Flag not handled");

		break;
	}

	#pragma omp parallel	
	Periodic_Constrain_Particles_To_Box();

	Gravity_Periodic_Init(); 

	Profile("Init");

	return ;
}

void Allocate_Particle_Structures()
{
	#pragma omp parallel // Task is threadprivate
	{

	const double npart_per_rank = (double) Sim.Npart_Total/(double) Sim.NRank;

	Task.Npart_Total_Max = npart_per_rank * PART_ALLOC_FACTOR;

	for (int i = 0; i < NPARTYPE; i++)
		Task.Npart_Max[i] = (double)Sim.Npart[i]/Sim.NRank * PART_ALLOC_FACTOR;

	} // omp parallel

	size_t nBytes = Task.Npart_Total_Max * sizeof(*P);

	P = Malloc(nBytes, "P");

	rprintf("\nReserving space for %jd particles per task in *P,"
			" factor %g\n", Task.Npart_Total_Max, PART_ALLOC_FACTOR);

	//G = Malloc(Task.Npart_Max[0] * sizeof(*G), "G");

#ifndef MEMORY_MANAGER
	memset(P, 0, Task.Npart_Total_Max * sizeof(*P) );
	//memset(G, 0, Task.Npart_Total_Max * sizeof(*G) );
#endif

	return ;
}

