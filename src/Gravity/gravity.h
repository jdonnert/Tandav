#ifdef GRAVITY

void Accel_Gravity_Simple();

#ifdef GRAVITY_TREE

struct Tree_Node {
	uint32_t Bitfield; 	// bit 0-5:level, 6-8:key, 9-31:free
	int DNext;		   	// distance to the next node; or -ipart
	Float Pos[3];		// Node Center
	Float CoM[3];		// Center of Mass
	Float Vel_CoM[3];		// Velocity of Center of Mass
	float Mass;			// Total Mass
	int DUp;			// distance to the parent
	int Npart;			// number of particles in node
} *Tree;

size_t NNodes;
size_t Max_Nodes;

void Gravity_Tree_Build();
void Gravity_Tree_Update();
void Gravity_Tree_Acceleration();

#endif // GRAVITY_TREE

#endif // GRAVITY
