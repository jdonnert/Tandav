#ifdef GRAVITY

#ifdef GRAVITY_SIMPLE
void Gravity_Simple_Accel();
#else
inline void Gravity_Simple_Accel(){};
#endif // GRAVITY_SIMPLE

#ifdef GRAVITY_TREE

extern struct Tree_Node {
	int DNext;			// Distance to the next node; or particle -DNext-1
	uint32_t Bitfield;	// bit 0-5:level, 6-8:key, 9:local, 10:top, 11-31:free
	int DUp;			// Number of nodes to the parent
	int Npart;			// Number of particles in node
	Float Pos[3];		// Node Center
	Float Mass;			// Total Mass of particles inside node
	Float CoM[3];		// Center of Mass
	Float Dp[3];		// Velocity of Center of Mass
} *Tree;

int NNodes, NTop_Nodes;

void Gravity_Tree_Build();
void Gravity_Tree_Acceleration();
void Gravity_Tree_Update_Kicks(const int ipart, const double dt);
void Gravity_Tree_Update_Topnode_Kicks();
void Gravity_Tree_Update_Drift(const double dt);

int Level(const int node); // bitfield functions

enum Tree_Bitfield { LOCAL=9, TOP=10, UPDATED=11 }; // offset by one

bool Node_Is(const enum Tree_Bitfield bit, const int node);
void Node_Set(const enum Tree_Bitfield bit, const int node);
void Node_Clear(const enum Tree_Bitfield bit, const int node);

#ifdef PERIODIC
void Gravity_Tree_Periodic();
void Gravity_Tree_Periodic_Init();
#else // !PERIODIC
inline void Gravity_Tree_Periodic(){};
inline void Gravity_Tree_Periodic_Init(){};
#endif // PERIODIC

#endif // GRAVITY_TREE

#ifdef GRAVITY_MULTI_GRID
void Gravity_Multi_Grid_Long_Range();
#else
inline void Gravity_Multi_Grid_Long_Range() {};
#endif // GRAVITY_MULTI_GRID

#else // ! GRAVITY

inline void Gravity_Tree_Update_Kicks(const int ipart, const double dt) {};
inline void Gravity_Tree_Update_Topnode_Kicks() {};
inline void Gravity_Tree_Update_Drift(const double dt) {};
inline void Gravity_Tree_Periodic(){};
inline void Gravity_Tree_Periodic_Init(){};
inline void Gravity_Multi_Grid_Long_Range() {};
inline void Accel_Gravity_Simple(){};

#endif // GRAVITY
