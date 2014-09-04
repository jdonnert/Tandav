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

/* This exposes the time integration of the code. Everything else is found
 * in Update() */

int main(int argc, char *argv[])
{
	preamble(argc, argv);

	Read_and_Init();

	Setup();
	
	Update(BEFORE_MAIN_LOOP);
		
	for (;;) { // run, Forest, run !
		
		Set_New_Timesteps();

		Kick_Halfstep();

		if (time_For_Snapshot()) {
		
			Drift_To_Snaptime();
			
			Write_Snapshot();
		}
 		
		Drift_To_Sync_Point();

		Make_Active_Particle_List();

		Update(BEFORE_FORCES);

		Compute_Forces();

		Kick_Halfstep();
	
		if (time_Is_Up())
			break;

	if (Sig.WriteRestartFile) 
		Write_Restart_File();

	}

	Finish();

	return EXIT_SUCCESS;
}

/* Here we do OpenMP and MPI init: 
 * Because of Thread Parallelism, every thread has a rank 
 * Task.Rank which is a combination of MPI rank and ThreadID. On 
 * every MPI rank there is a main thread on which
 * Task.IsThreadMain is true. Every thread is treated as its own
 * MPI task, including communication, still loops allow work-sharing */

static void preamble(int argc, char *argv[])
{
	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

	Assert(provided == MPI_THREAD_MULTIPLE, 
			"MPI thread multiple not supported, have %d", provided);

	MPI_Comm_rank(MPI_COMM_WORLD, &Task.MPI_Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Sim.NTask);

	MPI_Is_thread_main(&Task.IsThreadMain);

#pragma omp parallel
   	{
   		Task.ThreadID = omp_get_thread_num();
   		Sim.NThreads = omp_get_num_threads();
		
		Sim.NRank = Sim.NTask ; // * Sim.NThreads;
		
		Task.Rank = Task.MPI_Rank ;//* Sim.NThreads + Task.ThreadID;

		Task.Seed[2] = 1441981 * Task.ThreadID; // init thread safe std rng
	   	erand48(Task.Seed); // remove leading 0 in some implementations
	
		if (Task.MPI_Rank == MASTER && Task.ThreadID == MASTER)
			Task.IsMaster = true;
   	}

	if (Task.IsMaster) {

		printf("# Tandav #\n\n");

		printf("Using %d MPI tasks, %d OpenMP threads \n\n", 
				Sim.NTask, Sim.NThreads);
		
		Print_compile_time_settings();

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
		Param.StartFlag = atoi(argv[2]);

	if (Param.StartFlag == 10) 
		Write_Parameter_File(Param.File); // dead end

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}	

static bool time_Is_Up()
{
	if (Sig.Endrun) {
	
		rprintf("Encountered Signal: Endrun, t=%g", Time.Current);

		return true;
	}
	
	if (Time.IntCurrent == Time.IntEnd) {
		
		rprintf("EndTime reached: %g \n", Time.End);
		
		return true;
	}

	if (Runtime() >= Param.RuntimeLimit) {

		rprintf("Runtime limit reached: %g\n", Param.RuntimeLimit);

		Sig.WriteRestartFile = true;

		return true;
	}

	return false;
}

static bool time_For_Snapshot()
{
	if (Sig.WriteSnapshot) {

		Sig.WriteSnapshot = false;
	
		rprintf("Encountered Signal: Write Snapshot  t=%g \n", Time.Current);

		return true;
	}

	if (Time.Current + Time.Step >= Time.NextSnap) { 
	
		rprintf("Snapshot No. %d at t=%g, Next at t=%g \n", 
				Time.SnapCounter+1, Time.NextSnap, Time.NextSnap+Time.BetSnap);

		return true;
	}

	return false;
}
