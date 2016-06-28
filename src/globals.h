#ifndef GLOBALS_H
#define GLOBALS_H

/*
 * Unexposed Code Parameters
 */

#define CHARBUFSIZE 512L	// Maximum No. of bytes in every char buffer
#define NPARTYPE 6L			// No of particle types
#define MEM_ALIGNMENT 64L	// memory alignment
#define MASTER 0L			// Global master MPI task for printing

/*
 * Types
 */

#ifdef DOUBLE_PRECISION // your compiler must understand __uint128_t integers

typedef uint64_t intime_t;		// type of integer time 
typedef uint64_t ID_t;			// type of particle ID
typedef double Float;			// type of floating point variables
typedef __uint128_t peanoKey; 	// long peanokey, 128 bit = 42 triplets/levels
typedef uint64_t shortKey;		// short peanokey, 64 bit = 21 triplets/levels
#define MPI_MYFLOAT MPI_DOUBLE 	// corresponding MPI communication type macro
#define SQRT sqrt				// "type aware" square root function for speed
#define EXP exp					// "type aware" exponential function for speed

#else // ! DOUBLE_PRECISION    

typedef uint32_t intime_t;		// safe memory ...	
typedef uint32_t ID_t;		
typedef float Float;
typedef uint64_t peanoKey; 		
typedef uint32_t shortKey;		
#define MPI_MYFLOAT MPI_FLOAT
#define SQRT sqrtf				// tgmath.h is ugly ...
#define EXP expf

#endif // ! DOUBLE_PRECISION

extern void Finish();
extern void Print_Compile_Time_Settings();

/* 
 * Workaround
 */

double erand48(unsigned short xsubi[3]);

/*
 * Global variables
 */

extern int Master;				// Global Rank Master (for printing)
extern int NRank;				// Number of MPI tasks
extern int NThreads;			// Number of OpenMP threads
extern int NTask;				// NRank * NThreads

extern struct Local_Task_Properties {
	int ID;						// unique ID of thread
	int Rank;					// MPI Rank of this thread
	int Thread_ID;				// OpenMP ID of this thread
	bool Is_Master;				// == true on global master rank 
	bool Is_MPI_Master;			// == true on MPI master rank 
	bool Is_Thread_Main;		// == true on local thread masters
	int Npart[NPARTYPE];		// Number of particles on this rank
	int Npart_Total;			// Sum of Npart
	uint64_t Npart_Total_Max;	// per task taking into account imbalance.
	uint64_t Npart_Max[NPARTYPE];// Use this if array size scales with Npart
	size_t Buffer_Size;			// for Thread Safe Buffer
	unsigned short Seed[3];		// Thread safe urand48() seed
} Task;
#pragma omp threadprivate(Task) // modifications only in parallel env. !!

extern struct Global_Simulation_Properties {
	uint64_t Npart_Total;		// total global number of particles
	uint64_t Npart[NPARTYPE];	// global number of particles
	double Mpart[NPARTYPE];		// Global Masses from header
	double Boxsize[3];			// Now in 3D !
} Sim;

extern int * restrict Active_Particle_List;
int NActive_Particles;

enum Start_Parameters {
	READ_IC = 0,
	READ_RESTART = 1,
	READ_SNAP = 2,
	DUMP_PARFILE = 10
};
#endif // GLOBALS_H
