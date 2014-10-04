/* 
 * Initialise global variables 
 */

#include "globals.h"
#include "io/io.h"
#include "domain.h"

struct Parameters_From_File Param; 
struct Global_Simulation_Properties Sim;

struct Local_Task_Properties Task = { 0 };
#pragma omp threadprivate(Task)

struct Particle_Data *P = NULL;

static void init_particles();

void Read_and_Init() 
{
	Read_Parameter_File(Param.File);
	
 	Init_Profiler();

	Init_Memory_Management();

	Init_Logs();

	init_particles();

	Init_Domain_Decomposition();

	switch (Param.Start_Flag) {

	case 0: 

		Read_Snapshot(Param.Input_File);
		
		break;

	case 1: 
			
		Read_Restart_File();
			
		break;

	default: 
			
		Assert(0, "Start Flag not handled");

		break;
	}
	
	return ;
}

static void init_particles()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*P);

	P = Malloc(nBytes);

	return ;
}
