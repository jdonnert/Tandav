/* INCLUDES */
#include <stdlib.h> 		// system       
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <mpi.h> 		// parallelisation
#include <omp.h>

#include <gsl/gsl_math.h> 	// GNU scientific library
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

/* CODE PARAMETERS */
#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.3	// Mem overhead in P to account for inbalance

/* VARIABLES */
extern struct Local_Task_Properties {		
	int Rank;				// MPI Rank of this processor
	int ThreadID;			// OpenMP ID of this thread
	int Npart[NPARTYPE];	// Number of particles on this processor
	int NpartTotal;			// Sum of Npart
	unsigned short Seed[3];	// Thread local urand48() seed
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
	char File[CHARBUFSIZE]; 	// parameter file name
	char InputFile[CHARBUFSIZE];
	char OutputFileBase[CHARBUFSIZE];
	int StartFlag;				// invokation mode
	int NumIOTasks;				// written in parallel
	int MaxMemSize;				// Memory Ceiling in 1024^2 Bytes
	int NumOutputFiles;			// Number of output files
	float TimeLimit;			// Time Limit
} Param;

extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	float Mass;
	uint32_t ID;
	uint32_t TimeBin;
	uint64_t Peanokey;
} *P;
