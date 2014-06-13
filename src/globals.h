#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

/* CODE PARAMETERS */
#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.2	// Mem overhead to account for dynamic inbalance
#define MASTER 0			// Rank of MPI master task

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
	uint64_t NpartTotalMax;		// per task taking into account imbalance.
	uint64_t NpartMax[NPARTYPE];// Use this if array size scales with Npart
	double Mpart[NPARTYPE]; 	// Global Masses  from header
	double Boxsize[3];			// Now in 3D !
} Sim;

extern struct Parameters_From_File {
	char File[CHARBUFSIZE]; 	// parameter file name
	char InputFile[CHARBUFSIZE];
	char OutputFileBase[CHARBUFSIZE];
	int StartFlag;				// invokation mode
	int NumIOTasks;				// written in parallel
	int MaxMemSize;				// Memory Ceiling in 1024^2 Bytes
	int NumOutputFiles;			// Number of files per snapshot
	double RuntimeLimit;			// Runtime Limit
	double TimeIntAccuracy;		// Accuracy of Time integration
	double GravSoftening;		// Gravitational softening parameter
	int CommBufSize;			// in 1024 Bytes
} Param;

extern struct Simulation_Flags { // a primitive signal implementation
	bool Fullstep;				// Current step if Fullstep
	bool WriteSnapshot;			// write a snapshot
	bool WriteRestartFile;		// write a restart file upon exit
	bool Endrun;					// stops the runs regulary
} Flag;

extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	float Acc[3];
	uint32_t ID;
	int Type;
	float Mass;
	int TimeBin;
	peanoKey Peanokey;
#ifdef OUTPUT_FORCE
	float Force[3];
#endif // OUTPUT_FORCE
} *P;

#endif // GLOBALS_H
