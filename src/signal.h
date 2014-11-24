bool Time_Is_Up();
bool Time_For_Snapshot();
void Test_For_Domain_Update();

extern struct Simulation_Signals { // communicate an event across the code
	bool Fullstep;				// Current step is fullstep
	bool Write_Snapshot;		// write a snapshot this iteration
	bool Write_Restart_File;	// write a restart file upon exit
	bool Endrun;				// stops the run
	bool Drifted_To_Snaptime;	// results in integer timeline out of sync
	bool First_Step;			// doing first step of the simulation
	bool Domain_Update;			// do domain decomposition & tree build now
	bool Force_Domain;
	bool Force_Tree_Build;		
} Sig;
