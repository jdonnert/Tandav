#ifdef GRAVITY

void Accel_Gravity_Simple();

#ifdef GRAVITY_TREE

struct Tree_Node {
	uint32_t Bitfield; 	// bit 0-5:level, 6-8:key, 9-31:free
	int DNext;		   	// distance to the next node; or -ipart
	Float CoM[3];		// Center of Mass
	float Mass;			// Total Mass
	int DUp;			// distance to the parent
	int Npart;			// number of particles in node
#ifdef TREE_POSITIONS
	Float Pos[3];		// Node Center
#endif
} *Tree;

size_t NNodes;
size_t Max_Nodes;

void Build_Gravity_Tree();
void Gravity_Tree_Acceleration();

#endif // GRAVITY_TREE

#endif // GRAVITY
