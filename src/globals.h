#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

int * restrict Active_Particle_List, NActive_Particles;

extern struct Local_Task_Properties {		
	int Rank;				// MPI Rank of this thread
	int Thread_ID;			// OpenMP ID of this thread
	bool Is_Master;			// == true on global master rank 
	bool Is_MPI_Master;		// == true on MPI master rank 
	int Is_Thread_Main;		// == true on local thread masters
	int Npart_Total;		// Sum of Npart
	int Npart[NPARTYPE];	// Number of particles on this processor
	uint64_t Npart_Total_Max;	// per task taking into account imbalance.
	uint64_t Npart_Max[NPARTYPE];// Use this if array size scales with Npart
	unsigned short Seed[3];	// Thread safe urand48() seed
} Task;
#pragma omp threadprivate(Task) // modifications only in parallel env. !!

extern struct Global_Simulation_Properties {	
	int Master;					// Global Rank Master
	int NRank;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	int NTask;					// NRank * NThreads
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
	intime_t Int_Time_Pos;		// position on integer timeline
	float Cost;					// computational weight of particle
#ifdef GRAVITY_TREE
	int Tree_Parent;			// Tree node containing particle
#endif
	ID_t ID; // add below 
	Float Pos[3];
	Float Vel[3];
	Float Acc[3];
	Float Mass;
#ifdef GRAVITY_POTENTIAL
	Float Grav_Pot;
#endif
#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	Float Grav_Acc[3];
#endif
} *P;

#endif // GLOBALS_H
