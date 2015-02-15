
/*
 * The Domain decomposition creates bunches in *D that will later be 
 * transformed into top nodes. The tree walk during acceleration then 
 * traverses the topnode list and kicks off particle export 
 * or local tree walk.
 * The top nodes store a "Target" which is the node where the subtree starts in 
 * "*Tree", if it is positive. If "Target" is negative, it points to the MPI 
 * rank that holds the subtree. This corresponds to the "P[ipart].Tree_Parent" 
 * index, that point to the parent leave of particle ipart in *Tree if 
 * positive and to a top node index "-1*P[ipart].Tree_Parent - 1" if 
 * negative to avoid the degeneracy of index "0".
 */

union Domain_Node_List {

	struct Bunch_Node {	// Data needed for Domain Decomposition
		shortKey Key;		// Largest Peano key held by this bunch
		int Target;			// MPI rank
		int Level;
		size_t Npart;
		float Cost;			// cpu times
		int First_Part;		// starts the tree build
		bool Is_Local;		// on this rank
		int Modify;			// split  this bunch
	} Bunch;

#ifdef GRAVITY_TREE
	struct Top_Tree_Node {	//  dynamic top nodes, tree entry points
		shortKey Key;		// Number of nodes to the parent
		int Target;			// Tree/part index (>=0) or MPI rank - 1 (<0)
		int Level;			// Top node level
		int Npart;			// Number of particles in node
		float Pos[3];		// Node Center
		float Mass;			// Total Mass of particles inside node
		float CoM[3];		// Center of Mass
		float Dp[3];		// Velocity of Center of Mass, add above ! 
	} TNode;

#endif //GRAVITY_TREE

} *D;

int NBunches;

struct Domain_Properties { // smallest cubic box containing all particles
	double Size;
	double Origin[3];
	double Center[3];
	double Center_Of_Mass[3];
} Domain;

void Setup_Domain_Decomposition();
void Domain_Decomposition();
void Find_Global_Center_Of_Mass(double *CoM_out);
