#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

/* CODE PARAMETERS */
#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.2	// Mem overhead to account for dynamic inbalance

/* VARIABLES */
extern struct Local_Task_Properties {		
	int Rank;				// MPI Rank of this processor
	int ThreadID;			// OpenMP ID of this thread
	int NpartTotal;			// Sum of Npart
	int Npart[NPARTYPE];	// Number of particles on this processor
	unsigned short Seed[3];	// Thread safe urand48() seed
} Task;
#pragma omp threadprivate(Task)

extern struct Global_Simulation_Properties {	
	int NTask;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	uint64_t NpartTotal;		// total global number of particles
	uint64_t Npart[NPARTYPE]; 	// global number of particles
	uint64_t NpartTotalMean;	// As above, but taking into account imbalance.
	uint64_t NpartMean[NPARTYPE];// Use this if array size scales with Npart
	double Mpart[NPARTYPE]; 	// Global Masses  from header
	double Boxsize[3];			// Now in 3D !
	double CurrentTime;			// Holds current simulation time
} Sim;

extern struct Parameters_From_File {
	char File[CHARBUFSIZE]; 	// parameter file name
	char InputFile[CHARBUFSIZE];
	char OutputFileBase[CHARBUFSIZE];
	int StartFlag;				// invokation mode
	int NumIOTasks;				// written in parallel
	int MaxMemSize;				// Memory Ceiling in 1024^2 Bytes
	int NumOutputFiles;			// Number of files per snapshot
	float RuntimeLimit;			// Runtime Limit
	int CommBufSize;			// in 1024 Bytes
} Param;

extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	uint32_t ID;
	uint32_t TimeBin;
	peanoKey Peanokey;
	int Type;
	float Mass;
#ifdef OUTPUT_FORCE
	float Force[3];
#endif // OUTPUT_FORCE
} *P;

#endif // GLOBALS_H
