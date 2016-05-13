bool Time_Is_Up();
bool Time_For_Snapshot();
bool Time_For_Domain_Update();

extern struct Simulation_Signals { // communicate an event across the code
	bool Sync_Point;			// all particles synchronised
	bool Write_Snapshot;		// write a snapshot this iteration
	bool Write_Restart_File;	// write a restart file upon exit
	bool Endrun;				// stops the run
	bool Prepare_Step;			// preparing for the simulation
	bool First_Step;			// First step of the simulation
	bool Domain_Update;			// do domain decomposition & tree build now
	bool Tree_Update;			// use only with Domain_Update
	bool Use_BH_Criterion;		// Use different opening criterion
} Sig;
#pragma omp threadprivate(Sig)  // the compiler hates this to be public
