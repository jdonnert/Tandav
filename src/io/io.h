/* Parameter I/O */
struct parameter_definitions {
	char tag[CHARBUFSIZE]; // Parameter file tag
	char val[CHARBUFSIZE]; // Standard value
	void *addr; // Address of target variable
	enum param_type {
		FLOAT,
		INT,
		STRING,
	} type; // type of addr
};

static const struct parameter_definitions ParDef[] = { 
	{"Boxsize", "10000", &Param.Boxsize, FLOAT}, // tag, val, addr, type
	{"No_IOTasks", "1", &Param.No_Output_Files, INT},
	{"Input_File", "IC_file", &Param.Input_File, STRING},
	{"Output_File_Base", "snap_", &Param.Output_File_Base, STRING}
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
} ;

struct io_block_def {		// everything we need to define a Block in Format 2
	char label[4];
	char name[CHARBUFSIZE];
	enum target_variable {
		VAR_P, VAR_G
	} target;				// identify global var
	int offset;				// offset in underlying struct
	int nBytes; 			// over all elements of blockvar
	bool always_present; 	// can we live without it ?
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_MEMBER_SIZE(member) sizeof(((struct Particle_Data *)0)->member)

static const struct io_block_def Block[] = {
	{"POS ", "Positions", VAR_P, P_OFFSET(Pos), 3*P_MEMBER_SIZE(Pos[0]), 1},
	{"VEL ", "Velocities", VAR_P, P_OFFSET(Vel), 3*P_MEMBER_SIZE(Vel[0]), 1},
	{"ID  ", "IDs", VAR_P, P_OFFSET(ID), P_MEMBER_SIZE(ID), 1},
	{"MASS", "Masses", VAR_P, P_OFFSET(Mass), P_MEMBER_SIZE(Mass), 1}
};

#undef P_OFFSET
#undef P_MEMBER_SIZE

static const int NBlocks = sizeof(Block) / sizeof(*Block);

