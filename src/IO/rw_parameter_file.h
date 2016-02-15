
/* 
 * Parameter File I/O 
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
	{"BufferSize", "32", &Param.Buffer_Size, PAR_INT}, // in MB per MPI rank
	{"TimeLimitCPU", "20864", &Param.Runtime_Limit, PAR_DOUBLE},

	{"\n%% Simulation Characteristics %%\n", "", NULL, PAR_COMMENT},
	{"Boxsize0", "-1", &Sim.Boxsize[0], PAR_DOUBLE},
	{"Boxsize1", "-1", &Sim.Boxsize[1], PAR_DOUBLE},
	{"Boxsize2", "-1", &Sim.Boxsize[2], PAR_DOUBLE},
	{"TimeBegin", "0", &Time.Begin, PAR_DOUBLE},
	{"TimeEnd", "10", &Time.End, PAR_DOUBLE},
	{"TimeOfFirstSnaphot", "0", &Time.First_Snap, PAR_DOUBLE},
	{"TimeBetSnapshots", "0", &Time.Bet_Snap, PAR_DOUBLE},
	{"MaxSizeTimestep", "0.05", &Param.Max_Timestep, PAR_DOUBLE},
	{"MinSizeTimestep", "1e-7", &Param.Min_Timestep, PAR_DOUBLE},

	/* Add yours below */
};

static const int NTags = ARRAY_SIZE(ParDef);
