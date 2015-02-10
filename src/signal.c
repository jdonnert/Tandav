#include "globals.h"
#include "timestep.h"
#include "drift.h"

static bool test_for_stop_file();

#pragma omp threadprivate(Sig)
struct Simulation_Signals Sig;

/*
 * These functions handle simulation signals (->signal.h).
 */

bool Time_Is_Up()
{
	if (Sig.Endrun)
		rprintf("Encountered Signal: Endrun, t=%g", Time.Current);

	if (test_for_stop_file()) {
	
		rprintf("Found stop file t=%g\n", Time.Current);

		Sig.Write_Restart_File = true;

		Sig.Endrun = true;
	}
	
	if (Int_Time.Current == Int_Time.End) {

		rprintf("\nEndTime reached: %g \n", Time.End);
		
		Sig.Endrun = true;
	}

	if (Runtime() >= Param.Runtime_Limit) {

		rprintf("Runtime limit reached: t=%g at %g min\n", 
				Time.Current, Param.Runtime_Limit/60);

		Sig.Write_Restart_File = true;

		Sig.Endrun = true;
	}


	#pragma omp flush (Sig)

	return Sig.Endrun;
}

bool Time_For_Snapshot()
{
	if (Time.Current + Time.Step >= Time.Next_Snap) { 
		
		Drift_To_Snaptime();

		rprintf("\nSnapshot No. %d at t=%g, Next at t=%g \n", 
				Time.Snap_Counter+1, Time.Next_Snap, 
				Time.Next_Snap + Time.Bet_Snap);
	
		#pragma omp single
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

/*
 * Test if we have to do a domain decomposition & tree build, depending on
 * the number of interactions / drifted particles.
 */

static int Global_NPart_Updates = 0;
static int Local_NPart_Updates = 0;

void Check_For_Domain_Update()
{
	const double max_npart_updates = DOMAIN_UPDATE_PARAM*Sim.Npart_Total;

	Sig.Domain_Update = false;
	Sig.Tree_Update = false;

	#pragma omp single 
	{

	Local_NPart_Updates += NActive_Particles;
		
	MPI_Allreduce(&Local_NPart_Updates, &Global_NPart_Updates, 1, 
			MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	} // omp single

	#pragma omp flush (Global_NPart_Updates,Local_NPart_Updates)

	if (Sig.Fullstep || (Global_NPart_Updates > max_npart_updates)) {
		
		#pragma omp barrier

		#pragma omp single
		Global_NPart_Updates = Local_NPart_Updates = 0;
		
		Sig.Domain_Update = true;
		Sig.Tree_Update = true;
	}

	return ;
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

	//MPI_Bcast(&endrun, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	} // omp single

	#pragma omp flush

	return endrun;
}

