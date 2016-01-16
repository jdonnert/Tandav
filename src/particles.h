void Allocate_Particle_Structures();

struct Field_Def {	// define a particle property
	char Name[CHARBUFSIZE]; 	// 
	size_t Bytes; 	// sizeof member
	int N; 			// number of variables
};

#define P_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member)

static const struct Field_Def P_Fields[] = { 
	{"Type", 			P_SIZEOF(Type), 			1}
	,{"Time_Bin", 		P_SIZEOF(Time_Bin),			1}
	,{"Int_Time_Pos",	P_SIZEOF(Int_Time_Pos),		1}
	,{"ID", 		 	P_SIZEOF(ID),				1}
	,{"Cost",			P_SIZEOF(Cost),				1}
	,{"Pos", 			P_SIZEOF(Pos),				3}
	,{"Vel", 			P_SIZEOF(Vel),				3}
	,{"Acc", 			P_SIZEOF(Acc),				3}
	,{"Mass", 			P_SIZEOF(Mass),				1}
	,{"Grav_Acc",		P_SIZEOF(Grav_Acc),			3}
#ifdef GRAVITY_POTENTIAL
	,{"Grav_Pot",		P_SIZEOF(Grav_Pot),			1}
#endif	
#ifdef GRAVITY_TREE
	,{"Tree_Parent",	P_SIZEOF(Tree_Parent),		1}
#endif
	// Add yours here !
};

#undef P_FIELD_SIZE

extern const int NP_Fields;
