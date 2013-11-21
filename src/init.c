/* Initialise global variables */
#include "globals.h"
#include "proto.h"
#include "io/io.h"

struct Time_Integration_Infos Time;   
struct Particle_Data *P; 
struct Local_Task_Properties Task = { 0 };
struct Global_Simulation_Properties Sim;
struct Parameters_From_File Param; 
struct Unit_Constants_In_Cgs Unit = { 0 };

struct Particle_Data *P = NULL;

void Read_and_Init() 
{
	Read_Parameter_File(Param.File);
	
	Init_Memory_Management();

 	Init_Profiler();

	if (Param.StartFlag == 0) 
		Read_Snapshot(Param.InputFile);
	else if (Param.StartFlag == 1)
		Read_Restart_File();
	else 
		Assert(0, "Start Flag not handled");

#ifdef COMOVING
	Init_Cosmology();
#endif // COMOVING

	return ;
}

