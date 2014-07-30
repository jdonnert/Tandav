#include "globals.h"
#include "update.h"
#include "kick.h"
#include "drift.h"
#include "timestep.h"
#include "setup.h"
#include "peano.h"
#include "force.h"
#include "io/io.h"

static void preamble(int argc, char *argv[]);
static bool time_For_Snapshot();
static bool time_Is_Up();

/* This exposes the time integration */
int main(int argc, char *argv[])
{
	preamble(argc, argv);

	Read_and_Init();

	Setup();
			
	Update(BEFORE_MAIN_LOOP);
		
	Kick_Halfstep();
	
	for (;;) { // Quinn+97, Springel05
		
		Set_New_Timesteps();

		Kick_Halfstep();

		Drift();
		
		if (time_For_Snapshot())
			Write_Snapshot();

		if (time_Is_Up())
			break;

		Compute_Forces();

		Kick_Halfstep();
	}

	if (Flag.WriteRestartFile) 
		Write_Restart_File();

	Finish();

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
		
	Task.Seed[2] = 1441981 * Task.ThreadID; // init thread safe std rng
   	erand48(Task.Seed); // remove leading 0 in some implementations
   	}

	if (Task.Rank == MASTER) {

		printf("# Tandav #\n\n");

		printf("Using %d MPI tasks, %d OpenMP threads \n\n", 
				Sim.NTask, Sim.NThreads);
		
		Print_compile_time_settings();

		printf("\nsizeof(*P) = %zu byte\n", sizeof(*P));

		Assert(argc >= 2 && argc < 4, 
			"Wrong number of arguments, let me help you: \n\n" 
			"	USAGE: ./Tandav ParameterFile <StartFlag>\n\n"
			"	  0  : Read IC file and start simulation (default) \n"
			"	  1  : Read restart files and resume  \n"
			"	  2  : Read snapshot file and continue \n"
			"	  10 : Dump a valid parameter file for this Config\n");
	}

	strncpy(Param.File, argv[1], CHARBUFSIZE);
	
	if (argc > 2)
		Param.StartFlag = atoi(argv[2]);

	if (Param.StartFlag == 10) 
		Write_Parameter_File(Param.File); // dead end

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}	

static bool time_Is_Up()
{
	if (Time.IntCurrent == Time.IntEnd) {
		
		rprintf("EndTime reached: %g \n", Time.End);
		
		return true;
	}

	if (Flag.Endrun) 
		return true;

	if (Runtime() >= Param.RuntimeLimit) {

		rprintf("Runtime limit reached: %g\n", Param.RuntimeLimit);

		Flag.WriteRestartFile = true;

		return true;
	}

	return false;
}

static bool time_For_Snapshot()
{
	if (Flag.WriteSnapshot)
		return true;
 	
	if (Time.Current + Time.Step > Time.NextSnap) { // not exact
	
		Time.NextSnap += Time.BetSnap;

		rprintf("Snapshot at t=%g, Next Snapshot at t=%g \n", 
				Time.Current, Time.NextSnap);

		return true;
	}

	return false;
}
