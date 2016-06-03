#ifndef GRAVITY_TREE_H
#define GRAVITY_TREE_H

/*
 * Gravity using a conventional Barnes & Hutt tree algorithm
 */

#include "../includes.h"
#include "../domain.h"
#include "../periodic.h"

#if defined(GRAVITY) && defined(GRAVITY_TREE)

void Setup_Gravity_Tree();
void Gravity_Tree_Build();
void Gravity_Tree_Acceleration();
void Gravity_Tree_Update_Kicks(const int ipart, const double dt);
void Gravity_Tree_Update_Topnode_Kicks();
void Gravity_Tree_Update_Drift(const double dt);
void Gravity_Tree_Free();

extern struct Tree_Node {
	int DNext;			// Distance to the next node; or particle -DNext-1
	uint32_t Bitfield; 	// bit 0-5:level, 6-8:key, 9:local, 10:top
	int DUp;			// Number of nodes to the parent
	int Npart;			// Number of particles in node
	Float Pos[3];		// Node Center
	Float Mass;			// Total Mass of particles inside node
	Float CoM[3];		// Center of Mass
	Float Dp[3];		// Velocity of Center of Mass
} * restrict Tree;


double Epsilon[NPARTYPE], // softening
	   Epsilon2[NPARTYPE], 
	   Epsilon3[NPARTYPE];

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

#else // ! (GRAVITY && GRAVITY_TREE)

static inline void Setup_Gravity_Tree() {}; 
static inline void Gravity_Tree_Build() {};
static inline void Gravity_Tree_Acceleration() {};
static inline void Gravity_Tree_Periodic() {};
static inline void Gravity_Tree_Update_Kicks(const int ipart,const double dt){};
static inline void Gravity_Tree_Update_Topnode_Kicks() {};
static inline void Gravity_Tree_Update_Drift(const double dt) {};
static inline void Gravity_Tree_Free() {};

#endif // GRAVITY && GRAVITY_TREE

/* 
 * Periodic Boundaries need to walk the tree but interact with the
 * Ewald cube at the particle/node position.
 */

#if defined(GRAVITY) && defined(GRAVITY_TREE) && defined(PERIODIC)
void Gravity_Tree_Periodic();
void Tree_Periodic_Nearest(Float dr[3]);

#else

static inline void Gravity_Tree_Periodic() {};
static inline void Tree_Periodic_Nearest(Float dr[3]) {};
#endif // GRAVITY && GRAVITY_TREE && PERIODIC

#endif // GRAVITY_TREE_H
