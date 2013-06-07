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
#include "aux.h" 				// auxiliary functions (allocation, assert ...)
#include "cosmo.h" 				// cosmology functions 
#include "unit.h" 				// unit functions

/* CODE PARAMETERS */
#define CHARBUFSIZE 256 	// Maximum No. of char in every string
#define NO_PART_TYPES 6 		// No of particle types

/* VARIABLES */
extern struct Local_Task_Properties {		
	int Rank;					// MPI Rank of this processor
	int ThreadID;				// OpenMP ID of this thread
	int Npart[NO_PART_TYPES];	// Number of particles on this processor
	int NpartTotal;				// Sum of Npart
	int FirstActivePart;		// Start of linked list of active particles
	int *NextActivePart;		// Next in linked list
} Task;
#pragma omp threadprivate(Task)

extern struct Global_Simulation_Properties {	
	int NTask;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	uint64_t NpartTotal;		// total global number of particles
	uint64_t Npart[NO_PART_TYPES]; // global number of particles
	float Mpart[NO_PART_TYPES]; // Global Masses 
	float Boxsize;				
} Sim;

extern struct Parameters_From_File {
	int Start_Flag;				// invokation mode
	char File[CHARBUFSIZE]; 	// parameter file
	char Input_File[CHARBUFSIZE];
	char Output_File_Base[CHARBUFSIZE];
	int No_Output_Files;		// written in parallel
} Param;

extern struct Unit_Constants_In_Cgs {
	Length;
	Mass;
	Velocity;
	Time;
	Energy;
} Unit;

extern struct Time_Integration_Infos {
	float Begin;				// Start time of simulation
	float End;					// End time of simulation
	float Current;				// Current time of simulation
	float Base;					// Smallest step (Begin-End) / 2^63
	float NextSnap;				// Time of next snapshot to write
	bool Fullstep;				// Indicates if current step is fullstep
	size_t Nsteps;				// Number of steps walked so far
	int SnapCounter;			// Keep track of Snapshots written
} Time;

extern struct Particle_Data {
	float Pos[3];
	float Vel[3];
	uint64_t ID;
	float Mass;
} *P;

