void Allocate_Particle_Structures();

struct Field_Def {	// define a particle property
	size_t Offset; 	// address of member = address of P + Offset
	size_t Bytes; 	// sizeof member
	int N; 			// number of variables
};

#define P_OFFSET(member) offsetof(struct Particle_Data, member)
#define P_SIZEOF(member) sizeof(((struct Particle_Data *)0)->member)

static const struct Field_Def P_Fields[] = { 
	{P_OFFSET(Type), 			P_SIZEOF(Type), 			1}
	,{P_OFFSET(Time_Bin), 		P_SIZEOF(Time_Bin),			1}
	,{P_OFFSET(Int_Time_Pos),	P_SIZEOF(Int_Time_Pos),		1}
	,{P_OFFSET(ID), 		 	P_SIZEOF(ID),				1}
	,{P_OFFSET(Cost),			P_SIZEOF(Cost),				1}
	,{P_OFFSET(Pos), 			P_SIZEOF(Pos),				3}
	,{P_OFFSET(Vel), 			P_SIZEOF(Vel),				3}
	,{P_OFFSET(Acc), 			P_SIZEOF(Acc),				3}
	,{P_OFFSET(Mass), 			P_SIZEOF(Mass),				1}
	,{P_OFFSET(Grav_Acc),		P_SIZEOF(Grav_Acc),			3}
#ifdef GRAVITY_POTENTIAL
	,{P_OFFSET(Grav_Pot),		P_SIZEOF(Grav_Pot),			1}
#endif	
#ifdef GRAVITY_TREE
	,{P_OFFSET(Tree_Parent),		P_SIZEOF(Tree_Parent),		1}
#endif
	// Add yours here !
};

#undef P_OFFSET
#undef P_FIELD_SIZE

extern const int NP_Fields;
