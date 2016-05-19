#ifndef GLOBALS_H
#define GLOBALS_H

#include "includes.h" 

int Master;					// Global Rank Master
int NRank;					// Number of MPI tasks
int NThreads;				// Number of OpenMP threads
int NTask;					// NRank * NThreads

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
	double Total_Mass;			// sum over P.Mass, updated every timestep
	double Center_Of_Mass[3];	// sum of P.Mass*P.Pos
	double Kinetic_Energy;		// sum of 0.5 *P.Mass*P.Vel^2
	double Momentum[3];		    // sum of P.Mass*P.Vel
	double Angular_Momentum[3];	// sum of P.Mass*P.Pos x P.Vel
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

int * restrict Active_Particle_List;

extern struct Particle_Vector_Blocks{
	int * restrict First;
	int * restrict Last;
} V; // contingouos particle blocks on the same timestep

int NParticle_Vectors, NActive_Particles;

#endif // GLOBALS_H
