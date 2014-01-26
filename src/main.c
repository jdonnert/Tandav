#include "globals.h"
#include "proto.h"
#include "update.h"
#include "kick.h"
#include "drift.h"
#include "timestep.h"
#include "setup.h"
#include "io/io.h"

static void preamble(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	preamble(argc, argv);	

	Read_and_Init();

	Setup();

	Update(BEFORE_MAIN_LOOP);

	for (;;) {

		Kick_First_Halfstep();
		
		Update(AFTER_FIRST_KICK);

		Drift();

		if (Time.Current == Time.NextSnap)
			Write_Snapshot();

		if (Time.Current == Time.End
		 || Time.Running == Param.TimeLimit)
			break;

		Update(AFTER_DRIFT);

		Kick_Second_Halfstep();
		
		Update(AFTER_SECOND_KICK);
	}

	if (Time.Running == Param.TimeLimit) 
		Write_Restart_File();

	Finish();

	MPI_Finalize();

	return EXIT_SUCCESS;
}

static void preamble(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Task.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Sim.NTask);

#pragma omp parallel
   	{
   	Task.ThreadID = omp_get_thread_num();
   	Sim.NThreads = omp_get_num_threads();
		
	Task.Seed[2] = 1441981 * Task.ThreadID; // init thread safe rng
   	erand48(Task.Seed); // remove leading 0 in some implementations
   	}

	if (!Task.Rank) {
		printf("# Tandav #\n\n");

		printf("Using %d MPI tasks, %d OpenMP threads \n\n", 
				Sim.NTask, Sim.NThreads);
		
		Print_compile_time_settings();

		printf("\nsizeof(*P) = %zu byte\n", sizeof(*P));

		Assert(argc >= 2, 
			"Wrong number of arguments, let me help you: \n\n" 
			"USAGE: ./Tandav ParameterFile <StartFlag>\n\n"
			" 0  : Read IC file and start simulation (default) \n"
			" 1  : Read restart files and resume  \n"
			" 2  :  Read snapshot file and continue \n"
			" 10 : Dump a valid parameter file for this Config\n");
	}

	strncpy(Param.File, argv[1], CHARBUFSIZE);
	
	if (argc > 2)
		Param.StartFlag = atoi(argv[2]);

	if (Param.StartFlag == 10) 
		Write_Parameter_File(Param.File); // dead end

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}	
