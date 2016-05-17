void Allocate_Particle_Structures();
char * Select_Particle(const size_t field, const int comp, const int ipart);

/*
 * Here start the particle structures, which hold most of the data of the
 * code. Because we are using structures containing arrays, not an array of 
 * structures, automatic allocation etc need a description of these. These are
 * in P_Fields...
 */

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
	int * restrict Tree_Parent;	// Tree node leave, negative-1 if top node only
#endif
} P;

extern struct Gas_Particle_Data {
	Float * restrict Entropy;
	Float * restrict Volume;
	Float * restrict Density;
	Float * restrict Bfld[3];
} G;

extern struct Star_Particle_Data {
	Float * restrict Star_Formation_Rate;
} S;

extern struct Black_Hole_Particle_Data {
	Float * restrict Entropy;
} B;

extern size_t sizeof_P; // set in particles.c


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
	,{"Key",			P_SIZEOF(Key),				1} // keep first
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
