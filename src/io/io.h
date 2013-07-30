/* Parameter I/O */
struct parameter_definitions {
	char tag[CHARBUFSIZE]; // Parameter file tag
	char val[CHARBUFSIZE]; // Standard value
	void *addr; // Address of target variable
	enum param_type {
		FLOAT, INT, STRING,
	} type; // type of addr
};

static const struct parameter_definitions ParDef[] = { 
	{"InputFile", "IC_file", &Param.InputFile, STRING},
	{"OutputFileBase", "snap", &Param.OutputFileBase, STRING},
	{"NumIOTasks", "1", &Param.NumIOTasks, INT},
	{"NumOutputFiles", "1", &Param.NumOutputFiles, INT},
	{"Boxsize", "4000", &Sim.Boxsize, FLOAT},
	{"MaxMemSize", "1000", &Param.MaxMemSize, INT},
	{"Omega0", "1", &Cosmo.Omega0, FLOAT},
	{"OmegaLambda", "0.7", &Cosmo.OmegaLambda, FLOAT},
	{"HubbleParam", "0.7", &Cosmo.HubbleParam, FLOAT}
};

static const int NTags = sizeof(ParDef) / sizeof(*ParDef);

/* Snapshot I/O */
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
		VAR_P, VAR_G, VAR_S, VAR_BH, VAR_BND
	} Target;				// identify global var
	int Offset;				// offset in underlying struct
	int Nbytes; 			// sizeof target member
	int PartBitMask;		// == 1 at bit i+1, if needed by type i
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_MEMBER_SIZE(member) sizeof(((struct Particle_Data *)0)->member)

static const struct io_block_def Block[] = {
	{"POS ", "Positions", VAR_P, P_OFFSET(Pos), P_MEMBER_SIZE(Pos), 0x3F},
	{"VEL ", "Velocities", VAR_P, P_OFFSET(Vel), P_MEMBER_SIZE(Vel),0x3F},
	{"ID  ", "Short IDs", VAR_P, P_OFFSET(ID), P_MEMBER_SIZE(ID), 0x3F},
	{"MASS", "Masses", VAR_P, P_OFFSET(Mass), P_MEMBER_SIZE(Mass), 0x00}
};

#undef P_OFFSET
#undef P_MEMBER_SIZE

static const int NBlocks = sizeof(Block) / sizeof(*Block);

