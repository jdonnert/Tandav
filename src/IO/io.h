#include "../timestep.h"

/* 
 * Function Prototypes 
 */

extern void Read_Parameter_File(const char *);
extern void Write_Parameter_File(const char *);
extern void Read_Snapshot();
extern void Read_Restart_File();
extern void Read_and_Init();
extern void Write_Snapshot();
extern void Write_Restart_File();

/* Snapshot I/O */
unsigned int Largest_Block_Member_Nbytes();
unsigned int Npart_In_Block(const int, const int *);

struct gadget_header { // standard gadget header, filled to 256 byte
	uint32_t Npart[6];
	double Massarr[6];
	double Time;
	double Redshift;
	int32_t Flag_Sfr;
	int32_t Flag_Feedback;
	uint32_t Nall[6];
	int32_t Flag_Cooling;
	int32_t Num_Files;
	double Boxsize;
	double Omega0;
	double Omega_Lambda;
	double Hubble_Param;
	int32_t Flag_Age;
	int32_t Flag_Metals;
	uint32_t Nall_High_Word[6];
	char fill_bytes[59];
};

struct io_block_def {  // everything we need to define a Block in Format 2
	char Label[5];
	char Name[CHARBUFSIZE];
	enum target_variable {
		VAR_P,
		VAR_GAS,
		VAR_DM,
		VAR_STAR,
		VAR_DISK,
		VAR_BND
	} Target;			// identify global var
	size_t Offset;		// offset in underlying struct
	size_t Nbytes;		// sizeof target field
	bool IC_Required;	// needed on readin from ICs ?
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member)

static const struct io_block_def Block[] = {
	{"POS ", "Positions", VAR_P, P_OFFSET(Pos), P_SIZEOF(Pos), true},
	{"VEL ", "Velocities", VAR_P, P_OFFSET(Vel), P_SIZEOF(Vel),true},
	{"ID  ", "Short IDs", VAR_P, P_OFFSET(ID), P_SIZEOF(ID), true},
	{"MASS", "Masses", VAR_P, P_OFFSET(Mass), P_SIZEOF(Mass), false}

#ifdef OUTPUT_TOTAL_ACCELERATION
	,{"ACC", "Acceleration", VAR_P, P_OFFSET(Acc), P_SIZEOF(Acc), false}
#endif

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	,{"GACC", "Grav Acceleration", VAR_P, P_OFFSET(Grav_Acc),
		P_SIZEOF(Grav_Acc), false}
#endif

#ifdef OUTPUT_GRAV_POTENTIAL
	,{"GPOT", "Grav Potential", VAR_P, P_OFFSET(Grav_Pot), P_SIZEOF(Grav_Pot),
		false}
#endif

#ifdef OUTPUT_PEANO_KEY
	,{"PKEY","Peanokey",VAR_P,P_OFFSET(Peanokey),P_SIZEOF(peanoKey), false}
#endif

	// Add yours below 
};

#undef P_OFFSET
#undef P_FIELD_SIZE

static const int NBlocks = ARRAY_SIZE(Block);
