#ifndef DOMAIN_H
#define DOMAIN_H

/*
 * The Domain decomposition creates bunches in *D that are tree nodes at which
 * the domain will be decomposed. They will later be  transformed into "top 
 * nodes", below which are subtree will be build. They are actually useful for 
 * the tree walk, as they hold the particle CoM etc. The tree walk during 
 * acceleration then traverses the topnode list and kicks off particle export 
 * or local tree walk.
 * The top nodes store a "Target" which is the node where the subtree starts in 
 * "*Tree", if it is positive. If "Target" is negative, it points to the MPI 
 * rank that holds the subtree. This corresponds to the "P[ipart].Tree_Parent" 
 * index, that point to the parent leave of particle ipart in *Tree if 
 * positive and to a top node index "-1*P[ipart].Tree_Parent - 1" if 
 * negative to avoid the degeneracy of index "0".
 */


#include "includes.h"
#include "particles.h" 
#include "memory.h" 
#include "peano.h"
#include "sort.h"
#include "select.h"
#include "timestep.h"
#include "properties.h"

union Domain_Node_List {

	struct Bunch_Node {		// Data needed for Domain Decomposition
		shortKey Key;		// Largest Peano key held by this bunch
		int Target;			// MPI rank
		int Level;
		int First_Part;		// starts the tree build
		uint64_t Npart;
		float Cost;			// cpu times
		bool Is_Local;		// on this rank
		int Modify;			// split  this bunch
	} Bunch;

	struct Top_Tree_Node {	//  dynamic top nodes, tree entry points
		shortKey Key;		// Number of nodes to the parent
		int Target;			// Tree/part index (>=0) or MPI rank - 1 (<0)
		int Level;			// Top node level
		union {
			int First_Part;	// starts the tree build
			int First_Leaf; // save starting leave after tree build
		};
		union {
			int Npart;		// Number of particles in node (before tree_build)
			int Nleafs;		// Number of Leafs in node (after tree build)
		};
		float Pos[3];		// Node Center
		float Mass;			// Total Mass of particles inside node
#ifdef GRAVITY_TREE
		float CoM[3];		// Center of Mass
		float Dp[3];		// Velocity of Center of Mass, add above ! 
#endif //GRAVITY_TREE
	} TNode;

} * restrict D;

uint32_t NTop_Nodes;

struct Domain_Properties {
	double Size; // smallest cubic box containing all particles
	double Origin[3];
	double Center[3];
	double Center_Of_Mass[3];
} Domain;

#ifdef HIGHRES_REGION 
struct Domain_Properties Region[HIGHRES_REGION];
#endif

void Domain_Decomposition();
void Setup_Domain_Decomposition();
void Finish_Domain_Decomposition();

#endif // DOMAIN_H
