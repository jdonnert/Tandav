#include "globals.h"
#include "IO/io.h"
#include "Gravity/gravity.h"
#include "Gravity/gravity_periodic.h"

struct Parameters_From_File Param;
struct Global_Simulation_Properties Sim;

#pragma omp threadprivate(Task)
struct Local_Task_Properties Task = { 0 };

struct Particle_Data P = { NULL };
struct Gas_Particle_Data G = { NULL };

void Read_and_Init(int argc, char *argv[])
{
	Init_Profiler();

	Profile("Init");

	Read_Parameter_File(Param.File);

	Init_Memory_Management();

	Init_Logs();

	Init_Units();

	Init_Constants();

	Init_Cosmology(); // COMOVING

	/* Add yours above */

	switch (Param.Start_Flag) {

	case READ_IC:

		Read_Snapshot(Param.Input_File); // also init particle structures

		break;

	case READ_RESTART:

		Read_Restart_File();

		break;

	case READ_SNAP:

		Assert(argc > 3, "Missing snapshot number in program invokation");    

		Restart.Snap_Counter =  atoi(argv[3]);

		char snap_file[CHARBUFSIZE] = {""};

		sprintf(snap_file, "%s_%03d", Param.Output_File_Base, atoi(argv[3]));

		Read_Snapshot(snap_file);

		break;

	default:

		Assert(false, "Start Flag %d not handled", Param.Start_Flag);

		break;
	}

	Init_Periodic(); // PERIODIC

	#pragma omp parallel
	Periodic_Constrain_Particles_To_Box(); // PERIODIC
		
	Gravity_Periodic_Init();

	Profile("Init");

	return ;
}



