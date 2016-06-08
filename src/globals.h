#ifndef GLOBALS_H
#define GLOBALS_H

/*
 * Unexposed Code Parameters
 */

#define CHARBUFSIZE 512L	// Maximum No. of bytes in every char buffer
#define NPARTYPE 6L			// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define MASTER 0L			// Global master MPI task for printing

/*
 * Types
 */

#ifdef DOUBLE_PRECISION
typedef double Float;		// type of floating point variables in P
#define MPI_MYFLOAT MPI_DOUBLE // corresponding MPI communication type macro
#define SQRT sqrt			// type aware square root function for speed
#else // ! DOUBLE_PRECISION
typedef float Float;
#define MPI_MYFLOAT MPI_FLOAT
#define SQRT sqrtf
#endif // ! DOUBLE_PRECISION

#ifdef LONG_IDS
typedef uint64_t ID_t;		// type of particle ID
#else
typedef uint32_t ID_t;		
#endif // ! LONG_IDS

typedef uint32_t intime_t;		// type of integer time 

typedef uint64_t shortKey;		// short peanokey, 64 bit = 21 triplets/levels
typedef __uint128_t peanoKey; 	// long peanokey, 128 bit = 42 triplets/levels

enum Start_Parameters {
	READ_IC = 0,
	READ_RESTART = 1,
	READ_SNAP = 2,
	DUMP_PARFILE = 10
};

extern void Finish();
extern void Print_Compile_Time_Settings();

/* 
 * Workaround
 */

double erand48(unsigned short xsubi[3]);

/*
 * Global variables
 */

extern int Master;				// Global Rank Master
extern int NRank;				// Number of MPI tasks
extern int NThreads;			// Number of OpenMP threads
extern int NTask;				// NRank * NThreads

extern struct Local_Task_Properties {
	int ID;						// unique ID of thread
	int Rank;					// MPI Rank of this thread
	int Thread_ID;				// OpenMP ID of this thread
	bool Is_Master;				// == true on global master rank 
	bool Is_MPI_Master;			// == true on MPI master rank 
	int Is_Thread_Main;			// == true on local thread masters
	int Npart_Total;			// Sum of Npart
	int Npart[NPARTYPE];		// Number of particles on this processor
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

extern struct Parameters_From_File {
	char File[CHARBUFSIZE];		// parameter file name
	int Start_Flag;
	char Input_File[CHARBUFSIZE];
	char Output_File_Base[CHARBUFSIZE];
	char Log_File_Dir[CHARBUFSIZE];
	int Num_IO_Tasks;			// written in parallel
	int Max_Mem_Size;			// Memory Ceiling in 1024^2 Bytes
	int Buffer_Size;			// Total size of thread safe buffer in MB
	int Num_Output_Files;		// Number of files per snapshot
	double Runtime_Limit;		// in sec
	double Max_Timestep;		// largest timestep constraint
	double Min_Timestep;		// smallest timestep constraint
	double Part_Alloc_Factor;	// Allowed mem imbalance in Particles
	double Time_Int_Accuracy;	// 
	double Grav_Softening[NPARTYPE]; // gravitiational softening
} Param;

extern int * restrict Active_Particle_List;
int NActive_Particles;
double arr[10];
#endif // GLOBALS_H
