void Allocate_Particle_Structures();
char * Select_Particle(const size_t field, const int comp, const int ipart);

/*
 * P_Fields provides a description of P, so we can loop through the components
 * and allocate, add and remove particles.
 */

struct Field_Def {	// define a particle property
	char Name[CHARBUFSIZE]; 	// 
	size_t Bytes; 	// sizeof member
	int N; 			// number of variables
};

#define P_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member[0])
#define P_SIZEOF3(member) sizeof(((struct Particle_Data *)0)->member[0][0])

static const struct Field_Def P_Fields[] = { 
	{"Type", 			P_SIZEOF(Type), 			1}
	,{"Time_Bin", 		P_SIZEOF(Time_Bin),			1}
	,{"It_Drift_Pos",	P_SIZEOF(It_Drift_Pos),		1}
	,{"It_Kick_Pos",	P_SIZEOF(It_Kick_Pos),		1}
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
	,{"Last_Acc_Mag",	P_SIZEOF(Last_Acc_Mag),	1}
#endif
	// Add yours here !
};

#undef P_SIZEOF
#undef P_SIZEOF3

extern const int NP_Fields;
