#ifndef PARTICLES_H
#define PARTICLES_H

#include "includes.h"
#include "timestep.h"
#include "memory.h"
#include "IO/parameter_file.h"

void Allocate_Particle_Structures();
char * Select_Particle(const size_t field, const int comp, const int ipart);

/*
 * Here start the particle structures, which hold most of the data of the
 * code. Because we are using structures containing arrays, not an array of 
 * structures, automatic allocation needs a description of the structure. 
 * These are in P_Fields, which we use to loop through the members of P and
 * allocate, move etc ...
 */

struct Field_Def {			// define particle descriptor
	char Name[CHARBUFSIZE]; // id
	size_t Bytes; 			// sizeof member
	int N; 					// dimension
};

extern struct Particle_Data {
	int * restrict Type;				// keep first
	int * restrict Time_Bin;
	intime_t * restrict It_Drift_Pos;	// drift position on integer timeline
	intime_t * restrict It_Kick_Pos;	// kick position on integer timeline
	peanoKey * restrict Key;			// Reversed peano key
	ID_t * restrict ID; 					 
	Float * restrict Cost;				// computational weight of particle
	Float * restrict Pos[3];
	Float * restrict Vel[3];
	Float * restrict Acc[3];
	Float * restrict Mass;
	Float * restrict Grav_Acc[3];
	Float * restrict Last_Acc_Mag;		// Magnitude of Last Acc for tree force
#ifdef GRAVITY_POTENTIAL
	Float * restrict Grav_Pot;
#endif
#ifdef GRAVITY_TREE
	int * restrict Tree_Parent;			// Tree node leaf, negative if tnode
#endif
#ifdef GRAVITY_FMM
	int * restrict Leaf_Index;			// particle in Leaf2Nodes
#endif
} P;

const static struct Field_Def P_Fields[] = { 
	 {"Type", 			sizeof(int), 		1} // keep first
	,{"Time_Bin", 		sizeof(int),		1}
	,{"It_Drift_Pos",	sizeof(intime_t),	1}
	,{"It_Kick_Pos",	sizeof(intime_t),	1}
	,{"Key",			sizeof(peanoKey),	1} 
	,{"ID", 		 	sizeof(ID_t),		1}
	,{"Cost",			sizeof(Float),		1}
	,{"Pos", 			sizeof(Float),		3}
	,{"Vel", 			sizeof(Float),		3}
	,{"Acc", 			sizeof(Float),		3}
	,{"Mass", 			sizeof(Float),		1}
	,{"Grav_Acc",		sizeof(Float),		3}
	,{"Last_Acc_Mag",	sizeof(Float),		1}
#ifdef GRAVITY_POTENTIAL
	,{"Grav_Pot",		sizeof(Float),		1}
#endif	
#ifdef GRAVITY_TREE
	,{"Tree_Parent",	sizeof(int),		1}
#endif
#ifdef GRAVITY_FMM
	,{"Leaf_Index",		sizeof(int),		1}
#endif
	// Add yours here !
};


extern struct Gas_Particle_Data {
	Float * restrict Entropy;
	Float * restrict Volume;
	Float * restrict Density;
	Float * restrict Bfld[3];
	// Add yours here !
} G;

const static struct Field_Def G_Fields[] = {
	 {"Entropy",			sizeof(Float), 		1}
	,{"Volume",				sizeof(Float), 		1}
	,{"Density",			sizeof(Float),		1}
	,{"Bfld",				sizeof(Float),		3}
	// Add yours here !
};

extern struct Star_Particle_Data {
	Float * restrict Star_Formation_Rate;
	// Add yours here !
} S;

const static struct Field_Def S_Fields[] = {
	 {"Star_Formation_Rate",	sizeof(Float), 		1}
	// Add yours here !
};

extern struct Black_Hole_Particle_Data {
	Float * restrict Entropy;
	// Add yours here !
} BH;

const static struct Field_Def BH_Fields[] = {
	 {"Entropy",	sizeof(Float), 		1}
	// Add yours here !
};



size_t sizeof_P; 
const int NP_Fields;


#endif // PARTICLES_H
