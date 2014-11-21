#include "globals.h"
#include "timestep.h"
#include "drift.h"

static bool test_for_stop_file();

struct Simulation_Signals Sig;

/*
 * These functions handle simulation signals (->signal.h).
 */

bool Time_Is_Up()
{
	if (Sig.Endrun)
		rprintf("Encountered Signal: Endrun, t=%g", Time.Current);

	if (test_for_stop_file()) {
	
		rprintf("Found stop file t=%g", Time.Current);

		Sig.Write_Restart_File = true;

		Sig.Endrun = true;
	}
	
	if (Int_Time.Current == Int_Time.End) {
		
		rprintf("EndTime reached: %g \n", Time.End);
		
		Sig.Endrun = true;
	}

	if (Runtime() >= Param.Runtime_Limit) {

		rprintf("Runtime limit reached: t=%g at %g min\n", 
				Time.Current, Param.Runtime_Limit/60);

		Sig.Write_Restart_File = true;

		Sig.Endrun = true;
	}

	return Sig.Endrun;
}

bool Time_For_Snapshot()
{
	if (Time.Current + Time.Step >= Time.Next_Snap) { 
		
		Drift_To_Snaptime();

		rprintf("\nSnapshot No. %d at t=%g, Next at t=%g \n", 
				Time.Snap_Counter+1, Time.Next_Snap, 
				Time.Next_Snap + Time.Bet_Snap);
	
		#pragma omp single nowait
		Time.Next_Snap += Time.Bet_Snap;

		return true;
	}

	if (Sig.Write_Snapshot) {

		Sig.Write_Snapshot = false;
	
		rprintf("\nEncountered Signal: Write Snapshot %d at t=%g \n", 
				Time.Snap_Counter, Time.Current);

		return true;
	}


	return false;
}

bool Time_For_Domain_Update()
{
	Sig.Domain_Updated = false;

	if (1)
		return true;

	return false;
}

static int endrun = false;

static bool test_for_stop_file()
{
	#pragma omp single
	{

	if (Task.Is_MPI_Master) {
			
		FILE *fp = fopen("./stop", "r");
			
		if (fp != NULL) {
			
			fclose(fp);

			endrun = true;
		}
	}

	} // omp single

	int result = 0;

	#pragma omp single
	MPI_Allreduce(&endrun, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	return (bool) result;
}

