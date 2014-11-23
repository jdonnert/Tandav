bool Time_Is_Up();
bool Time_For_Snapshot();
bool Time_For_Domain_Update();

extern struct Simulation_Signals { // communicate an event across the code
	bool Fullstep;				// Current step is fullstep
	bool Write_Snapshot;		// write a snapshot this iteration
	bool Write_Restart_File;	// write a restart file upon exit
	bool Endrun;				// stops the run
	bool Drifted_To_Snaptime;	// results in integer timeline out of sync
	bool First_Step;			// doing first step of the simulation
	bool Domain_Updated;		// did domain decomposition this step
	bool Force_Domain;
	bool Force_Tree_Build;		
} Sig;
