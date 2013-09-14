/* Intialise global variables */
#include "globals.h"
#include "proto.h"

struct Time_Integration_Infos Time;   
struct Particle_Data *P; 
struct Local_Task_Properties Task = { 0 };
struct Global_Simulation_Properties Sim;
struct Parameters_From_File Param; 
struct Unit_Constants_In_Cgs Unit = { 0 };

struct Particle_Data *P = NULL;

void Init() 
{
	Init_Memory_Management();

 	Init_Profiler();

	return ;
}

