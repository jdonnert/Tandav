/* 
 * Initialise global variables 
 */

#include "globals.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "domain.h"

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
		
		Sig.First_Step = Sig.Fullstep = true;

		break;

	case 1: 
			
		Read_Restart_File();
			
		break;

	default: 
			
		Assert(0, "Start Flag not handled");

		break;
	}
	
#ifdef PERIODIC
	Constrain_Particles_To_Box();
#endif
	
	Print_Memory_Usage();
	
	Profile("Init");

	return ;
}

void Allocate_Particle_Structures()
{
	const double npart_per_rank = (double)Sim.Npart_Total/(double) Sim.NRank;

	#pragma omp parallel // Task is threadprivate
	{
	
	Task.Npart_Total_Max = ceil(npart_per_rank * PARTALLOCFACTOR);

	for (int i = 0; i < NPARTYPE; i++)
		Task.Npart_Max[i] = ceil((double)Sim.Npart[i] / (double)Sim.NRank 
				* PARTALLOCFACTOR);
	
	} // omp parallel

	rprintf("\nReserving space for %llu particles per task, factor %g\n", 
			Task.Npart_Total_Max, PARTALLOCFACTOR);

	P = Malloc(Task.Npart_Total_Max * sizeof(*P), "P"); 

	//G = Malloc(Task.Npart_Max[0] * sizeof(*G), "G");

	return ;
}

