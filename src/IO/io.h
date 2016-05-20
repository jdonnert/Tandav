#ifndef IO_H
#define IO_H

#include "../includes.h"
#include "../cosmology.h"
#include "../timestep.h"

void Read_Parameter_File(const char *);
void Write_Parameter_File(const char *);
void Read_Snapshot();
void Read_Restart_File();
void Read_and_Init();
void Write_Snapshot();
void Write_Restart_File();



/* 
 * Snapshot I/O 
 */

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

/*
 * Block provides a description of all output blocks that can be written.
 * This way we only have to edit one place to add a block.
 */

struct io_block_def {  // everything we need to define a Block in Format 2
	char Label[5];
	enum target_variable {
		VAR_P,
		VAR_GAS,
		VAR_DM,
		VAR_STAR,
		VAR_DISK,
		VAR_BND
	} Target;			// identify global var
	size_t Offset;		// offset in underlying struct
	size_t Ncomp;		// vector length / number of components
	size_t Nbytes;		// sizeof target field
	bool IC_Required;	// needed on readin from ICs ?
	char Name[CHARBUFSIZE];
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)

const static struct io_block_def Block[] = {

	 {"POS ", VAR_P, P_OFFSET(Pos ), 3, sizeof(Float), true, "Positions"}
	,{"VEL ", VAR_P, P_OFFSET(Vel ), 3, sizeof(Float), true, "Velocities"}
	,{"ID  ", VAR_P, P_OFFSET(ID  ), 1, sizeof(ID_t), true, "Short IDs"}
	,{"MASS", VAR_P, P_OFFSET(Mass), 1, sizeof(Float), false, "Masses"}

#ifdef OUTPUT_TOTAL_ACCELERATION
	,{"ACC ", VAR_P, P_OFFSET(Acc ), 1, sizeof(Float), false, "Acceleration"}
#endif
#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	,{"GACC", VAR_P, P_OFFSET(Grav_Acc), 3, sizeof(Float), false, "Grav Accel"}
#endif
#ifdef OUTPUT_GRAV_POTENTIAL
	,{"GPOT", VAR_P, P_OFFSET(Grav_Pot), 1, sizeof(Float), false, "Grav Pot"}
#endif
#ifdef OUTPUT_PEANO_KEY
	,{"PKEY", VAR_P, P_OFFSET(Peanokey), 1, sizeof(PeanoKey), false, "Pkeys"}
#endif

	// Add yours below 
};

static const int NBlocks = ARRAY_SIZE(Block);

#endif // IO_H
