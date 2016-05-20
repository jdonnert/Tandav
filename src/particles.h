#ifndef PARTICLES_H
#define PARTICLES_H

#include "includes.h"

void Allocate_Particle_Structures();
char * Select_Particle(const size_t field, const int comp, const int ipart);

struct Field_Def {	
	char Name[CHARBUFSIZE]; // id
	size_t Bytes; 			// sizeof member
	int N; 					// dimension
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member)
#define P_SIZEOF3(member) sizeof(((struct Particle_Data *)0)->member[0][0])

const static struct Field_Def P_Fields[] = { 
	{"Type", 			P_SIZEOF(Type), 			1} // keep first
	,{"Time_Bin", 		P_SIZEOF(Time_Bin),			1}
	,{"It_Drift_Pos",	P_SIZEOF(It_Drift_Pos),		1}
	,{"It_Kick_Pos",	P_SIZEOF(It_Kick_Pos),		1}
	,{"Key",			P_SIZEOF(Key),				1} 
	,{"ID", 		 	P_SIZEOF(ID),				1}
	,{"Cost",			P_SIZEOF(Cost),				1}
	,{"Pos", 			P_SIZEOF3(Pos),				3}
	,{"Vel", 			P_SIZEOF3(Vel),				3}
	,{"Acc", 			P_SIZEOF3(Acc),				3}
	,{"Mass", 			P_SIZEOF(Mass),				1}
	,{"Grav_Acc",		P_SIZEOF3(Grav_Acc),		3}
#ifdef GRAVITY_POTENTIAL
	,{"Grav_Pot",		P_SIZEOF(Grav_Pot),			1}
#endif	
#ifdef GRAVITY_TREE
	,{"Tree_Parent",	P_SIZEOF(Tree_Parent),		1}
	,{"Last_Acc_Mag",	P_SIZEOF(Last_Acc_Mag),		1}
#endif
	// Add yours here !
};

size_t sizeof_P; 
const int NP_Fields;


#endif // PARTICLES_H
