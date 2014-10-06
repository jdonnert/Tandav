bool Time_Is_Up();
bool Time_For_Snapshot();

extern struct Simulation_Signals { // communicate an event across the code
	bool Fullstep;				// Current step is fullstep
	bool Write_Snapshot;		// write a snapshot this iteration
	bool Write_Restart_File;	// write a restart file upon exit
	bool Endrun;				// stops the run
	bool Drifted_To_Snaptime;
} Sig;
#pragma omp threadprivate(Sig)
