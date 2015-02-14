#ifdef GRAVITY

void Accel_Gravity_Simple();

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

int NNodes;

void Gravity_Tree_Build();
void Gravity_Tree_Acceleration();
void Gravity_Tree_Update_Kicks(const int ipart, const double dt);
void Gravity_Tree_Update_Topnode_Kicks();
void Gravity_Tree_Update_Drift(const double dt);
void Gravity_Tree_Ewald_Correction();
void Gravity_Tree_Short_Range();

int Level(const int node); // bitfield functions

enum Tree_Bitfield { LOCAL=9, TOP=10, UPDATED=11 }; // offset by one

bool Node_Is(const enum Tree_Bitfield bit, const int node);
void Node_Set(const enum Tree_Bitfield bit, const int node);
void Node_Clear(const enum Tree_Bitfield bit, const int node);

#endif // GRAVITY_TREE

#ifdef GRAVITY_MULTI_GRID

void Gravity_Multi_Long_Range();

#endif // GRAVITY_MULTI_GRID

#else // ! GRAVITY

void Gravity_Tree_Update_Kicks(const int ipart, const double dt) {};
void Gravity_Tree_Update_Topnode_Kicks() {};
void Gravity_Tree_Update_Drift(const double dt) {};

#endif // GRAVITY
