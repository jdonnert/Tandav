/*
 * Gravity using a conventional tree algorithm
 */

#if defined(GRAVITY) && defined(GRAVITY_TREE)

uint32_t NNodes;

struct Walk_Data_Particle { // stores exported particle data
	ID_t ID;
	Float Pos[3];
	Float Acc; 				// only magnitude of the last acceleration
	Float Mass;
};

struct Walk_Data_Result { 	// stores exported summation results
	Float Cost;
	double Grav_Acc[3];
#ifdef GRAVITY_POTENTIAL
	double Grav_Potential;
#endif
};

int Level(const int node); // bitfield functions

enum Tree_Bitfield { LOCAL=9, TOP=10, UPDATED=11 }; // offset by one

Float Node_Size(const int node);
bool Node_Is(const enum Tree_Bitfield bit, const int node);
void Node_Set(const enum Tree_Bitfield bit, const int node);
void Node_Clear(const enum Tree_Bitfield bit, const int node);

#endif // GRAVITY && GRAVITY_TREE

#if defined(GRAVITY) && defined(GRAVITY_TREE) && defined(PERIODIC)
void Gravity_Tree_Periodic();
void Tree_Periodic_Nearest(Float dr[3]);

#else

inline void Gravity_Tree_Periodic() {};
inline void Tree_Periodic_Nearest(Float dr[3]) {};
#endif // GRAVITY && GRAVITY_TREE && PERIODIC
