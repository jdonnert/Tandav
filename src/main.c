#include "globals.h"
#include "update.h"
#include "kick.h"
#include "drift.h"
#include "timestep.h"
#include "setup.h"
#include "domain.h"
#include "peano.h"
#include "accel.h"
#include "IO/io.h"

static void preamble(int argc, char *argv[]);

/* 
 * This exposes the time integration of the code. Everything else is found
 * in Update(). 
 */

int main(int argc, char *argv[])
{
	preamble(argc, argv);
	
	Read_and_Init();

	Setup();
	
	#pragma omp parallel
	{
	
	Update(BEFORE_MAIN_LOOP);

	#pragma omp barrier

	for (;;) { // run, Forest, run !
break;
		if (Time_Is_Up())
			break;

		Update(BEFORE_STEP);

		Set_New_Timesteps();

		Kick_First_Halfstep();

		if (Time_For_Snapshot()) 
			Write_Snapshot();
 		
		Drift_To_Sync_Point();
		
		Domain_Decomposition();

		Compute_Acceleration();

		Kick_Second_Halfstep();
	}

	} // omp parallel 

	if (Sig.Write_Restart_File) 
		Write_Restart_File();

	Finish();

	return EXIT_SUCCESS;
}

/* 
 * Here we do OpenMP and MPI init: 
 * Because of Thread Parallelism, every thread has a MPI rank 
 * Task.Rank and a ThreadID.
 * There is a global MPI master with Task.Is_MPI_Master = true. 
 * On every MPI rank there is a main thread on which 
 * Task.Is_Thread_Main = true. Every thread 
 * is treated as its own MPI task, including communication.
 * Always use Task.MPI_Rank inside an omp single region.
 */

static void preamble(int argc, char *argv[])
{
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	Assert(provided == MPI_THREAD_MULTIPLE, 
			"MPI thread multiple not supported, have %d", provided);

	MPI_Comm_rank(MPI_COMM_WORLD, &Task.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Sim.NRank);

	MPI_Is_thread_main(&Task.Is_Thread_Main);

	#pragma omp parallel
   	{
   	
		Task.Thread_ID = omp_get_thread_num();
   		Sim.NThreads = omp_get_num_threads();
		
		Sim.NTask = Sim.NRank * Sim.NThreads; 

		if (Task.Rank == MASTER && Task.Thread_ID == MASTER)
			Task.Is_Master = true;
		
		if (Task.Rank == MASTER)
			Task.Is_MPI_Master = true;
		
		Task.Seed[2] = 1441981L * Task.Thread_ID; // init thread safe std rng
	   	erand48(Task.Seed); // remove first 0 in some implementations

	} // omp parallel

	if (Task.Is_Master) {

		printf("# Tandav #\n\n");

		printf("Using %d MPI tasks, %d OpenMP threads \n", 
				Sim.NRank, Sim.NThreads);
		
		printf("\nsizeof(*P) = %zu byte\n", sizeof(*P)*CHAR_BIT/8);

		Assert(argc >= 2 && argc < 4, 
			"Wrong number of arguments, let me help you: \n\n" 
			"	USAGE: ./Tandav ParameterFile <StartFlag>\n\n"
			"	  0  : Read IC file and start simulation (default) \n"
			"	  1  : Read restart files and resume  \n"
			"	  2  : Read snapshot file and continue \n"
			"	 10  : Dump a valid parameter file for this Config\n");
	}

	strncpy(Param.File, argv[1], CHARBUFSIZE);
	
	if (argc > 2) // Start Flag given
		Param.Start_Flag = atoi(argv[2]);

	if (Param.Start_Flag == 10) 
		Write_Parameter_File(Param.File); // dead end

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}	


