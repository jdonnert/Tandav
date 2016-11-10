#include "signal.h"

static bool test_for_stop_file();

#pragma omp threadprivate(Sig)
struct Simulation_Signals Sig;

/* These functions handle simulation signals (->signal.h). */

bool Time_Is_Up()
{
	if (Sig.Endrun)
		rprintf("\nEncountered Signal: Endrun, t=%g\n\n", Time.Current);

	if (test_for_stop_file()) {

		rprintf("\nFound stop file t=%g\n\n", Time.Current);

		Sig.Restart_Write_File = true;

		Sig.Endrun = true;
	}

	if (Int_Time.Current == Int_Time.End) {

		rprintf("\nEndTime reached: %g \n\n", Time.End);

		Sig.Endrun = true;
	}

	return Sig.Endrun;
}

bool Runtime_Limit_Reached()
{
	if (Runtime() >= Param.Runtime_Limit) {

		rprintf("\nRuntime limit reached: t=%g at %g min\n\n",
				Time.Current, Param.Runtime_Limit/60);

		Sig.Restart_Write_File = true;

		Sig.Endrun = true;
	}

	return Sig.Endrun;
}

bool Time_For_Snapshot()
{
	if (Time.Current + Time.Step >= Time.Next_Snap) {
		
		rprintf("\nSnapshot No. %d at t=%g, Next at t=%g \n",
				Time.Snap_Counter, Time.Next_Snap,
				Time.Next_Snap + Time.Bet_Snap);

		#pragma omp barrier // no barrier in rprintf

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

/* Test if we have to do a domain decomposition & tree build, depending on
 * the number of interactions / drifted particles. */

static int Global_NPart_Updates = 0;
static int Local_NPart_Updates = 0;

bool Time_For_Domain_Update()
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

	if (Sig.Sync_Point || (Global_NPart_Updates > max_npart_updates)) {

		#pragma omp barrier

		#pragma omp single
		Global_NPart_Updates = Local_NPart_Updates = 0;

		Sig.Domain_Update = true;
		Sig.Tree_Update = true;
	}

	return Sig.Domain_Update;
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

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
