#include "../timestep.h"

/* Function Prototypes */
extern void Read_Parameter_File(const char *);
extern void Write_Parameter_File(const char *);
extern void Read_Snapshot();
extern void Read_Restart_File();
extern void Read_and_Init();
extern void Write_Snapshot();
extern void Write_Restart_File();

/* Parameter I/O */
typedef struct {
	char tag[CHARBUFSIZE]; // Parameter file tag
	char val[CHARBUFSIZE]; // Standard value
	void *addr; // Address of target variable
	enum param_type {
		DOUBLE, 
		INT, 
		STRING,
		COMMENT,
	} type; // type of addr
}  parameter;

static const parameter ParDef[] = { 

	{"\n%% Files %%\n", "", NULL, COMMENT},
	{"InputFile", "IC_file", &Param.InputFile, STRING},
	{"OutputFileBase", "snap", &Param.OutputFileBase, STRING},
	{"NumIOTasks", "1", &Param.NumIOTasks, INT},
	{"NumOutputFiles", "1", &Param.NumOutputFiles, INT},

	{"\n%% Code Parameters %%\n", "", NULL, COMMENT},
	{"MaxMemSize", "1000", &Param.MaxMemSize, INT},
	{"TimeLimitCPU", "20864", &Param.RuntimeLimit, INT},
	{"CommBufSize", "1000", &Param.CommBufSize, INT},

	{"\n%% Simulation Characteristics %%\n", "", NULL, COMMENT},
	{"Boxsize0", "4000", &Sim.Boxsize[0], DOUBLE},
	{"Boxsize1", "4000", &Sim.Boxsize[1], DOUBLE},
	{"Boxsize2", "4000", &Sim.Boxsize[2], DOUBLE},
	{"TimeBegin", "0", &Time.Begin, DOUBLE},
	{"TimeEnd", "10", &Time.End, DOUBLE},
	{"TimeOfFirstSnaphot", "0", &Time.FirstSnap, DOUBLE},
	{"TimeBetSnapshots", "0", &Time.BetSnap, DOUBLE},

	{"\n%% Cosmology %%\n", "", NULL, COMMENT},
	{"Omega0", "1", &Cosmo.Omega0, DOUBLE},
	{"OmegaLambda", "0.7", &Cosmo.OmegaLambda, DOUBLE},
	{"OmegaBaryon", "1", &Cosmo.OmegaBaryon, DOUBLE},
	{"HubbleParam", "0.7", &Cosmo.HubbleParam, DOUBLE},

	/* Add yours below */
};

static const int NTags = ARRAY_SIZE(ParDef);

/* Snapshot I/O */
unsigned int largest_block_member_nbytes();
unsigned int npart_in_block(const int, const int *);

struct gadget_header { // standard gadget header, filled to 256 byte
	uint32_t Npart[6];
	double Massarr[6];
	double Time;
	double Redshift;
	int32_t FlagSfr;
	int32_t FlagFeedback;
	uint32_t Nall[6];
	int32_t FlagCooling;
	int32_t NumFiles;
	double Boxsize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int32_t FlagAge;
	int32_t FlagMetals;
	uint32_t NallHighWord[6];
	char fill_bytes[59];
};

struct io_block_def {  // everything we need to define a Block in Format 2
	char Label[5];
	char Name[CHARBUFSIZE];
	enum target_variable {
		VAR_P, 
		VAR_G, 
		VAR_S, 
		VAR_BH, 
		VAR_BND
	} Target;		// identify global var
	size_t Offset;		// offset in underlying struct
	size_t Nbytes; 	// sizeof target field
	int PartBitMask;// == 1 at bit i+1, if required for type i
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_FIELD_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member)

static const struct io_block_def Block[] = {

  	{"POS ", "Positions", VAR_P, P_OFFSET(Pos), P_FIELD_SIZEOF(Pos), 0xFF},
  	{"VEL ", "Velocities", VAR_P, P_OFFSET(Vel), P_FIELD_SIZEOF(Vel),0xFF},
  	{"ID  ", "Short IDs", VAR_P, P_OFFSET(ID), P_FIELD_SIZEOF(ID), 0xFF},
  	{"MASS", "Masses", VAR_P, P_OFFSET(Mass), P_FIELD_SIZEOF(Mass), 0xFF}
#ifdef OUTPUT_FORCE
  	,{"FRCE", "Forces", VAR_P, P_OFFSET(Force), P_FIELD_SIZEOF(Force), 0xFF}
#endif
#ifdef OUTPUT_PEANO_KEY
  	,{"PKEY","Peanokey",VAR_P,P_OFFSET(Peanokey),P_FIELD_SIZEOF(peanoKey),0xFF}
#endif

	/* Add yours below */
};

#undef P_OFFSET
#undef P_FIELD_SIZE

static const int NBlocks = ARRAY_SIZE(Block);
