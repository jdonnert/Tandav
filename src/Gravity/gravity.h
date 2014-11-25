#ifdef GRAVITY

void Accel_Gravity_Simple();

#ifdef GRAVITY_TREE

struct Tree_Node {
	uint32_t Bitfield; 	// bit 0-5:level, 6-8:key, 9:local, 10:top, 11-31:free
	int DNext;		   	// Distance to the next node; or particle -DNext-1
	int DUp;			// Number of nodes to the parent
	Float Pos[3];		// Node Center
	Float CoM[3];		// Center of Mass
	Float Dp[3];		// Velocity of Center of Mass
	float Mass;			// Total Mass of particles inside node
	int Npart;			// Number of particles in node
} *Tree;

int NNodes;
int Max_Nodes;

void Gravity_Tree_Build();
void Gravity_Tree_Acceleration();
void Gravity_Tree_Update_Kicks(const Float* dp, const int node);
void Gravity_Tree_Update_Drift(const double dt);
void Gravity_Tree_Update_Topnode_Kicks();
void Gravity_Tree_Ewald_Correction();
void Gravity_Tree_Short_Range();

int Level(const int node); // bitfield functions

enum Tree_Bitfield { LOCAL=10, TOP, UPDATED };
bool Node_Is(const enum Tree_Bitfield bit, const int node);
void Node_Set(const enum Tree_Bitfield bit, const int node);
void Node_Clear(const enum Tree_Bitfield bit, const int node);

#endif // GRAVITY_TREE

#endif // GRAVITY
