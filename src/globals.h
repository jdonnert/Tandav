#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

typedef float Float;

/* CODE PARAMETERS */

#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.2	// Mem overhead for dynamic inbalance
#define MASTER 0			// Global master MPI thread

/* VARIABLES */

extern struct Local_Task_Properties {		
	bool IsMaster;			// == true on global master rank 
	int IsThreadMain;		// == true on local thread masters
	int Rank;				// combined OMP & MPI Rank of this processor
	int MPI_Rank;			// MPI Rank of this processor
	int ThreadID;			// OpenMP ID of this thread
	int NpartTotal;			// Sum of Npart
	int Npart[NPARTYPE];	// Number of particles on this processor
	unsigned short Seed[3];	// Thread safe urand48() seed
} Task;
#pragma omp threadprivate(Task)

extern struct Global_Simulation_Properties {	
	int Master;					// MPI Master Rank
	int NRank;					// NTask * NThreads
	int NTask;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	uint64_t NpartTotal;		// total global number of particles
	uint64_t Npart[NPARTYPE]; 	// global number of particles
	uint64_t NpartTotalMax;		// per task taking into account imbalance.
	uint64_t NpartMax[NPARTYPE];// Use this if array size scales with Npart
	double Mpart[NPARTYPE]; 	// Global Masses  from header
	double Boxsize[3];			// Now in 3D !
} Sim;

int * restrict ActiveParticleList, NActiveParticles;

extern struct Parameters_From_File {
	char File[CHARBUFSIZE]; 	// parameter file name
	char InputFile[CHARBUFSIZE];
	char OutputFileBase[CHARBUFSIZE];
	int StartFlag;				// invokation mode
	int NumIOTasks;				// written in parallel
	int MaxMemSize;				// Memory Ceiling in 1024^2 Bytes
	int NumOutputFiles;			// Number of files per snapshot
	double RuntimeLimit;		
	int CommBufSize;			// in 1024 Bytes
} Param;

extern struct Simulation_Signals { // communicate an event across the code
	bool Fullstep;				// Current step is fullstep
	bool WriteSnapshot;			// write a snapshot this iteration
	bool WriteRestartFile;		// write a restart file upon exit
	bool Endrun;				// stops the runs regularly
	bool Synchronize_Drift;		// drift to next sync point on integer timeline
} Sig;

extern struct Particle_Data {
	Float Pos[3];
	Float Vel[3];
	Float Force[3];
	Float Potential;
	Float Mass;
	uint32_t ID;
	int Type;
	uint32_t TimeBin;
	peanoKey Peanokey;
	/* add below */
} *P;

#endif // GLOBALS_H
