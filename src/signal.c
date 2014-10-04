#include "globals.h"
#include "timestep.h"
#include "drift.h"

#pragma omp threadprivate(Sig)
struct Simulation_Signals Sig;

bool Time_Is_Up()
{
	if (Sig.Endrun) {
	
		rprintf("Encountered Signal: Endrun, t=%g", Time.Current);

		return true;
	}
	
	if (Int_Time.Current == Int_Time.End) {
		
		rprintf("EndTime reached: %g \n", Time.End);
		
		return true;
	}

	if (Runtime() >= Param.Runtime_Limit) {

		rprintf("Runtime limit reached: %g min\n", Param.Runtime_Limit/60);

		Sig.Write_Restart_File = true;

		return true;
	}

	return false;
}

bool Time_For_Snapshot()
{
	if (Sig.Write_Snapshot) {

		Sig.Write_Snapshot = false;
	
		rprintf("Encountered Signal: Write Snapshot %d at t=%g \n", 
				Time.Snap_Counter, Time.Current);

		return true;
	}

	if (Time.Current + Time.Step >= Time.Next_Snap) { 
	
		Drift_To_Snaptime();

		rprintf("Snapshot No. %d at t=%g, Next at t=%g \n", 
				Time.Snap_Counter+1, Time.Next_Snap, 
				Time.Next_Snap + Time.Bet_Snap);
	
		#pragma omp single nowait
		Time.Next_Snap += Time.Bet_Snap;

		return true;
	}

	return false;
}
