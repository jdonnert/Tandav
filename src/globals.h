#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

typedef float Float;		// type of floating point variables in P
typedef uint32_t ID_t;		// type of particle ID
typedef uint64_t intime_t; 	// type of integer time 

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
	int Master;					// Global Rank Master
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
	char Log_File_Dir[CHARBUFSIZE];
	int Start_Flag;				// invokation mode
	int Num_IO_Tasks;			// written in parallel
	int Max_Mem_Size;			// Memory Ceiling in 1024^2 Bytes
	int Num_Output_Files;		// Number of files per snapshot
	int Comm_Buf_Size;			// in 1024 Bytes
	double Runtime_Limit;		// in sec
} Param;

extern struct Particle_Data {
	int Type;
	int Time_Bin;
	intime_t Int_Time_Pos;
	peanoKey Peanokey;
	Float Pos[3];
	Float Vel[3];
	Float Acc[3];
	Float Mass;
#ifdef GRAVITY_POTENTIAL
	Float Grav_Pot;
#endif
	ID_t ID; // add below 
} *P;

#endif // GLOBALS_H
