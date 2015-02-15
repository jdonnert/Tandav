#include "globals.h"
#include "update.h"
#include "kick.h"
#include "drift.h"
#include "timestep.h"
#include "setup.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "peano.h"
#include "accel.h"
#include "IO/io.h"

static void preamble(int argc, char *argv[]);

/* 
 * This exposes the time integration of the code. 
 * We use the HOLD integrator from Pelupessy+ 2012. 
 */

int main(int argc, char *argv[])
{
	preamble(argc, argv);

	Read_and_Init();

	Setup();

	#pragma omp parallel
	{

	Update(BEFORE_MAIN_LOOP);

	for (;;) { // run, Forest, run !

		if (Time_Is_Up())
			break;

		Update(BEFORE_STEP);

		Set_New_Timesteps();

		Update(BEFORE_FIRST_KICK);

		Kick_First_Halfstep();

		if (Time_For_Snapshot())
			Write_Snapshot();

		Drift_To_Sync_Point();

		if (Time_For_Domain_Update())
			Domain_Decomposition();

		Update(BEFORE_FORCES);

		Compute_Acceleration();

		Kick_Second_Halfstep();

		Update(AFTER_STEP);
	}

	exit_trap:;

	if (Sig.Write_Restart_File)
		Write_Restart_File();
	else
		Write_Snapshot();

	} // omp parallel 

	Finish();

	return EXIT_SUCCESS;
}

/* 
 * Here we do OpenMP and MPI init: 
 * Because of Thread Parallelism, every thread has an MPI rank 
 * Task.Rank and a Thread ID Task.Thread_ID.
 * There is a global MPI master with Task.Is_MPI_Master == true. 
 * On every MPI rank there is a main thread on which 
 * Task.Is_Thread_Main == true. Every thread 
 * is treated as its own MPI task, including communication.
 * Always use Task.MPI_Rank inside an omp single region. In a parallel region
 * use Task.Rank to identify a thread on an CPU.
 */

static void preamble(int argc, char *argv[])
{
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	Assert(provided == MPI_THREAD_MULTIPLE,
			"MPI thread multiple not supported, have %d :-(",
			provided);

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

		Task.Seed[2] = 14041981L * (Task.Thread_ID); // init thread safe rng

		erand48(Task.Seed); // remove first 0 in some implementations

	} // omp parallel

	if (Task.Is_Master) {

		printf("#### Tandav ####\n\n");

		Assert( (argc >= 2) && (argc < 4),
			"Wrong number of arguments, let me help you: \n\n"
			"	USAGE: ./Tandav ParameterFile <StartFlag>\n\n"
			"	  0  : Read IC file and start simulation (default) \n"
			"	  1  : Read restart files and resume  \n"
			"	  2  : Read snapshot file and continue \n"
			"	 10  : Dump a valid parameter file for this Config\n");

		Print_compile_time_settings();

		printf("\nsizeof(*P) = %zu byte\n", sizeof(*P)*CHAR_BIT/8);
		printf("sizeof(*D) = %zu byte\n", sizeof(*D)*CHAR_BIT/8);
#ifdef GRAVITY_TREE
		printf("sizeof(*Tree) = %zu byte\n", sizeof(*Tree)*CHAR_BIT/8);
#endif

		printf("\nUsing %d MPI tasks, %d OpenMP threads \n\n",
				Sim.NRank, Sim.NThreads);

	}

	strncpy(Param.File, argv[1], CHARBUFSIZE);

	if (argc > 2) // Start Flag given ?
		Param.Start_Flag = atoi(argv[2]);
	else
		Param.Start_Flag = 0;

	//if (Param.Start_Flag == 1) 
		//Restore_From_Restart_File();

	if (Param.Start_Flag == 10) {

		Write_Parameter_File(Param.File); // dead end

		Finish();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return ;

}


