#include "init.h"

struct Global_Simulation_Properties Sim;
int * restrict Active_Particle_List;

#pragma omp threadprivate(Task)
struct Local_Task_Properties Task = { 0 };

int Master = 0, NRank = 0, NThreads = 0, NTask = 0;

/*
 * Initialization before IC reading
 */

void Init(int argc, char *argv[])
{
	Init_Profiler();

	Profile("Init");

	Read_Parameter_File(Param.File);

	Init_Memory_Management();

	Init_Logs();

	Init_Units();

	Init_Constants();

	Init_Cosmology(); // COMOVING

	Profile("Init");

	return ;
}



