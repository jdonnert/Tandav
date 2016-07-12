#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include "../includes.h"
#include "../particles.h" 

/* 
 * Parameter File definition 
 */

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
	double Part_Alloc_Factor;	// Allowed mem imbalance in Particles
	double Grav_Softening[NPARTYPE]; // gravitiational softening
	double Time_Begin;
	double Time_End;
	double Time_First_Snap;
	double Time_Bet_Snap;
	double Time_Int_Accuracy;
	double Time_Max_Timestep;	// largest timestep constraint
	double Time_Min_Timestep;	// smallest timestep constraint
} Param;

typedef struct {
	char tag[CHARBUFSIZE]; // Parameter file tag
	char val[CHARBUFSIZE]; // Standard value
	void *addr; // Address of target variable
	enum param_type {
		PAR_DOUBLE,
		PAR_INT,
		PAR_STRING,
		PAR_COMMENT,
	} type; // type of addr
} parameter_def;

static const parameter_def ParDef[] = {

	{"\n%% Files %%\n", "", NULL, PAR_COMMENT},
	{"InputFile", "IC_file", &Param.Input_File, PAR_STRING},
	{"OutputFileBase", "snap", &Param.Output_File_Base, PAR_STRING},
	{"LogFileDir", "./log", &Param.Log_File_Dir, PAR_STRING},
	{"NumIOTasks", "1", &Param.Num_IO_Tasks, PAR_INT},
	{"NumOutputFiles", "1", &Param.Num_Output_Files, PAR_INT},

	{"\n%% Code Parameters %%\n", "", NULL, PAR_COMMENT},
	{"MaxMemSize", "1024", &Param.Max_Mem_Size, PAR_INT},
	{"BufferSize", "32", &Param.Buffer_Size, PAR_INT}, // in MB all threads
	{"PartAllocFactor", "1.1", &Param.Part_Alloc_Factor, PAR_DOUBLE},
	{"TimeLimitCPU", "20864", &Param.Runtime_Limit, PAR_DOUBLE},

	{"\n%% Simulation Characteristics %%\n", "", NULL, PAR_COMMENT},
#ifdef PERIODIC_NO_CUBE
	{"Boxsize_0", "-1", &Sim.Boxsize[0], PAR_DOUBLE}, 
	{"Boxsize_1", "-1", &Sim.Boxsize[1], PAR_DOUBLE}, 
	{"Boxsize_2", "-1", &Sim.Boxsize[2], PAR_DOUBLE}, 
#else
	{"Boxsize", "-1", &Sim.Boxsize[0], PAR_DOUBLE}, // -1 = get from snapshot
#endif
	{"TimeBegin", "0", &Param.Time_Begin, PAR_DOUBLE},
	{"TimeEnd", "10", &Param.Time_End, PAR_DOUBLE},
	{"TimeOfFirstSnaphot", "0", &Param.Time_First_Snap, PAR_DOUBLE},
	{"TimeBetSnapshots", "0", &Param.Time_Bet_Snap, PAR_DOUBLE},
	{"TimeIntAccuracy", "0.2", &Param.Time_Int_Accuracy, PAR_DOUBLE},
	{"MaxSizeTimestep", "0.05", &Param.Time_Max_Timestep, PAR_DOUBLE},
	{"MinSizeTimestep", "1e-7", &Param.Time_Min_Timestep, PAR_DOUBLE},
#ifdef GRAVITY
	{"\n%% Gravity %%\n", "", NULL, PAR_COMMENT},
	{"GravSoftening", "10", &Param.Grav_Softening[1], PAR_DOUBLE}, // Plummer
#endif
	/* Add yours below */
};

static const int NTags = ARRAY_SIZE(ParDef);

#endif // PARAMETER_FILE_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
