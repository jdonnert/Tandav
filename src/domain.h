
/*
 * The Domain decomposition creates bunches in *D that will later be made
 * into the top nodes. The walk then traverses the topnode list and kicks
 * off particle export or local tree walk.
 */

union Domain_Node_List { 
	
	struct Bunch_Node { // Data needed for Domain Decomposition
		shortKey Key;	// Largest Peano key held by this bunch
		int Target; 	// MPIRANK
		int Npart;
		int Level;
		float Cost;
		int First_Part;
	} Bunch;

#ifdef GRAVITY_TREE
	struct Top_Tree_Node {	//  dynamic top nodes, tree entry points
		shortKey Key;		// Number of nodes to the parent
		int Target;	   		// Tree/part index (>0) or MPI rank (<0)
		int Npart;			// Number of particles in node
		float Pos[3];		// Node Center
		float Mass;			// Total Mass of particles inside node
		float CoM[3];		// Center of Mass
		float Dp[3];		// Velocity of Center of Mass, add above ! 
	} TNode; 
#endif //GRAVITY_TREE

} *D;


struct Domain_Properties { // smallest cubic box containing all particles
	double Size;	 
	double Origin[3];	
	double Center[3];	
} Domain;

int NBunches;
int NLocal_Bunches;

void Domain_Decomposition();
void Init_Domain_Decomposition();



