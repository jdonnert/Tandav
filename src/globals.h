/* INCLUDES */
#include <stdlib.h> 			// system       
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <mpi.h> 				// parallelisation
#include <omp.h>

#include <gsl/gsl_math.h> 		// GNU scientific library
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

#include "config.h" 			// holds configuration switches

#include "macro.h" 				// macro definitions
#include "aux.h" 				// auxiliary functions 
#include "cosmo.h" 				// cosmology functions 
#include "unit.h" 				// unit functions
#include "memory.h"				// memory management

/* CODE PARAMETERS */
#define CHARBUFSIZE 256L 		// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 			// No of particle types
#define MEM_ALIGNMENT 128L		// byte memory alignment
#define PARTALLOCFACTOR 1.3		// Mem overhead in P to account for inbalance

/* VARIABLES */
extern struct Local_Task_Properties {		
	int Rank;					// MPI Rank of this processor
	int ThreadID;				// OpenMP ID of this thread
	int Npart[NPARTYPE];		// Number of particles on this processor
	int NpartTotal;				// Sum of Npart
} Task;
#pragma omp threadprivate(Task)

extern struct Global_Simulation_Properties {	
	int NTask;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	uint64_t NpartTotal;		// total global number of particles
	uint64_t Npart[NPARTYPE]; 	// global number of particles
	float Mpart[NPARTYPE]; 		// Global Masses 
	float Boxsize;	
} Sim;

extern struct Parameters_From_File {
	int StartFlag;				// invokation mode
	char File[CHARBUFSIZE]; 	// parameter file name
	char InputFile[CHARBUFSIZE];
	char OutputFileBase[CHARBUFSIZE];
	int NumIOTasks;				// written in parallel
	int MaxMemSize;				// Memory Ceiling in 1024^2 Bytes
	int NumOutputFiles;			// Number of output files
} Param;

extern struct Time_Integration_Infos {
	float Begin;				// Start time of simulation
	float End;					// End time of simulation
	float Current;				// Current time of simulation
	float Base;					// Smallest step (Begin-End) / 2^63
	float NextSnap;				// Time of next snapshot to write
	bool Fullstep;				// Indicates if current step is fullstep
	int Nsteps;					// Number of steps walked so far
	int SnapCounter;			// Keep track of Snapshots written
	float Running;				// Run time of this task
	float Limit;				// Time Limit
} Time;

extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	uint32_t ID;
	float Mass;
} *P;
