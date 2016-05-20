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

extern struct Particle_Vector_Blocks{
	int * restrict First;
	int * restrict Last;
} V; // contingouos particle blocks on the same timestep

int NParticle_Vectors, NActive_Particles;

/*
 * Here start the particle structures, which hold most of the data of the
 * code. Because we are using structures containing arrays, not an array of 
 * structures, automatic allocation needs a description of the structure. 
 * These are in P_Fields, which we use to loop through the members of P and
 * allocate, move etc ...
 */

extern struct Particle_Data {
	int * restrict Type;				// keep first
	int * restrict Time_Bin;
	intime_t * restrict It_Drift_Pos;	// drift position on integer timeline
	intime_t * restrict It_Kick_Pos;	// kick position on integer timeline
	peanoKey * restrict Key;			// Reversed peano key
	ID_t * restrict ID; 					 
	Float * restrict Cost;				// computational weight of particle
	Float * restrict Pos[3];
	Float * restrict Vel[3];
	Float * restrict Acc[3];
	Float * restrict Mass;
	Float * restrict Grav_Acc[3];
	Float * restrict Last_Acc_Mag;		// Magnitude of Last Acc for tree force
#ifdef GRAVITY_POTENTIAL
	Float * restrict Grav_Pot;
#endif
#ifdef GRAVITY_TREE
	int * restrict Tree_Parent;			// Tree node leave, negative-1 if
#endif									// top node only

} P;


extern struct Gas_Particle_Data {
	Float * restrict Entropy;
	Float * restrict Volume;
	Float * restrict Density;
	Float * restrict Bfld[3];
} G;

extern struct Star_Particle_Data {
	Float * restrict Star_Formation_Rate;
} S;

extern struct Black_Hole_Particle_Data {
	Float * restrict Entropy;
} B;

#endif // GLOBALS_H
