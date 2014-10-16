/* 
 * Initialise global variables 
 */

#include "globals.h"
#include "io/io.h"
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
		
	Print_compile_time_settings();

	Read_Parameter_File(Param.File);

	Init_Memory_Management();

	Init_Logs();

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
	
#ifdef PERIODIC
	Constrain_Particles_To_Box();
#endif

	Profile("Init");

	return ;
}
