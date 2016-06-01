#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include "../includes.h"
#include "../timestep.h"

/* 
 * Parameter File definition 
 */

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
}  parameter;

static const parameter ParDef[] = {

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
	{"GravSoftening", "10", &Param.Grav_Softening[1], PAR_DOUBLE}, // Plummer

	{"\n%% Simulation Characteristics %%\n", "", NULL, PAR_COMMENT},
#ifdef PERIODIC_NO_CUBE
	{"Boxsize_0", "-1", &Sim.Boxsize[0], PAR_DOUBLE}, 
	{"Boxsize_1", "-1", &Sim.Boxsize[1], PAR_DOUBLE}, 
	{"Boxsize_2", "-1", &Sim.Boxsize[2], PAR_DOUBLE}, 
#else
	{"Boxsize", "-1", &Sim.Boxsize[0], PAR_DOUBLE}, 
#endif
	{"TimeBegin", "0", &Time.Begin, PAR_DOUBLE},
	{"TimeEnd", "10", &Time.End, PAR_DOUBLE},
	{"TimeOfFirstSnaphot", "0", &Time.First_Snap, PAR_DOUBLE},
	{"TimeBetSnapshots", "0", &Time.Bet_Snap, PAR_DOUBLE},
	{"TimeIntAccuracy", "0.2", &Param.Time_Int_Accuracy, PAR_DOUBLE},
	{"MaxSizeTimestep", "0.05", &Param.Max_Timestep, PAR_DOUBLE},
	{"MinSizeTimestep", "1e-7", &Param.Min_Timestep, PAR_DOUBLE},

	/* Add yours below */
};

static const int NTags = ARRAY_SIZE(ParDef);

#endif // PARAMETER_FILE_H
