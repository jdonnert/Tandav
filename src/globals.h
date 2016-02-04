#ifndef GLOBALS_H
#define GLOBALS_H

#include "proto.h"

int * restrict Active_Particle_List, NActive_Particles;

extern struct Local_Task_Properties {
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
	int Master;					// Global Rank Master
	int NRank;					// Number of MPI tasks
	int NThreads;				// Number of OpenMP threads
	int NTask;					// NRank * NThreads
	uint64_t Npart_Total;		// total global number of particles
	uint64_t Npart[NPARTYPE];	// global number of particles
	double Mpart[NPARTYPE];		// Global Masses from header
	double Boxsize[3];			// Now in 3D !
	double Total_Mass;			// sum over P.Mass, updated every timestep
	double Center_Of_Mass[3];	// center of Mass, updated every timestep
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
} Param;

/*
 * Here start the particle structures, which hold most of the data of the
 * code. Because we are using structures containing arrays, not an array of 
 * structures, automatic allocation etc need a description of these. You can
 * find these descriptors in particles_fields.h, the allocator in particles.c
 * We might auto generate these later.
 */

extern struct Particle_Data {
	int * restrict Type;				// keep this first !
	int * restrict Time_Bin;
	intime_t * restrict Int_Time_Pos;	// current position on integer timeline
	ID_t * restrict ID; 					 
	Float * restrict Cost;				// computational weight of particle
	Float * restrict Pos[3];
	Float * restrict Vel[3];
	Float * restrict Acc[3];
	Float * restrict Mass;
	Float * restrict Grav_Acc[3];
#ifdef GRAVITY_POTENTIAL
	Float * restrict Grav_Pot;
#endif
#ifdef GRAVITY_TREE
	int * restrict Tree_Parent;	// Tree node leave, negative-1 if top node only
#endif
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

extern size_t sizeof_P; // set in particles.c

#endif // GLOBALS_H
