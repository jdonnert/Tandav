#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

typedef float Float;
typedef uint32_t IDtype;

/* CODE PARAMETERS */

#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.2	// Mem overhead for dynamic inbalance
#define MASTER 0			// Global master MPI thread

/* VARIABLES */

int * restrict Active_Particle_List, NActive_Particles;

extern struct Local_Task_Properties {		
	bool Is_Master;			// == true on global master rank 
	int Is_Thread_Main;		// == true on local thread masters
	int Rank;				// combined OMP & MPI Rank of this thread
	int MPI_Rank;			// MPI Rank of this thread
	int Thread_ID;			// OpenMP ID of this thread
	int Npart_Total;		// Sum of Npart
	int Npart[NPARTYPE];	// Number of particles on this processor
	uint64_t Npart_Total_Max;	// per task taking into account imbalance.
	uint64_t Npart_Max[NPARTYPE];// Use this if array size scales with Npart
	unsigned short Seed[3];	// Thread safe urand48() seed
} Task;
#pragma omp threadprivate(Task)

extern struct Global_Simulation_Properties {	
	int Master;					// MPI Master Rank
	int NRank;					// NTask * NThreads
	int NTask;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	uint64_t Npart_Total;		// total global number of particles
	uint64_t Npart[NPARTYPE]; 	// global number of particles
	double Mpart[NPARTYPE]; 	// Global Masses  from header
	double Boxsize[3];			// Now in 3D !
} Sim;

extern struct Parameters_From_File {
	char File[CHARBUFSIZE]; 	// parameter file name
	char Input_File[CHARBUFSIZE];
	char Output_File_Base[CHARBUFSIZE];
	int Start_Flag;				// invokation mode
	int Num_IO_Tasks;			// written in parallel
	int Max_Mem_Size;			// Memory Ceiling in 1024^2 Bytes
	int Num_Output_Files;		// Number of files per snapshot
	double Runtime_Limit;		
	int Comm_Buf_Size;			// in 1024 Bytes
} Param;

extern struct Simulation_Signals { // communicate an event across the code
	bool Fullstep;				// Current step is fullstep
	bool Write_Snapshot;		// write a snapshot this iteration
	bool Write_Restart_File;	// write a restart file upon exit
	bool Endrun;				// stops the runs regularly
	bool Synchronize_Drift;		// drift to next sync point on integer timeline
} Sig;

extern struct Particle_Data {
	int Type;
	int Time_Bin;
	peanoKey Peanokey;
	Float Pos[3];
	Float Vel[3];
	Float Acc[3];
	Float Potential;
	Float Mass;
	IDtype ID;
	/* add below */
} *P;

#endif // GLOBALS_H
