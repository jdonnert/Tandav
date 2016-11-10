#ifndef FMM_H
#define FMM_H

#include "../includes.h"
#include "../domain.h"

/* Gravity using the Fast Multipole Method with dual tree traversal.
 * We use the notation from Dehnen & Read 2011 */

#ifdef GRAVITY_FMM

#define FMM_EXP_ORDER 1
#define NM FMM_EXP_ORDER 		// number of multipole components
#define NFN 3*(FMM_EXP_ORDER) // number of field tensor components

void Gravity_FMM_Acceleration();
void Gravity_FMM_Setup();
void Gravity_FMM_Free();
void Gravity_FMM_Build_P2M();
void Gravity_FMM_M2L();
void Gravity_FMM_P2L();
void Gravity_FMM_Update_Kicks();
void Gravity_FMM_Update_Topnode_Kicks();
void Gravity_FMM_Update_Drift(const double dt);

void M2L_Setup();

struct FMM_Node { // Dehnen & Read 2012

	uint32_t * restrict Bitfield;	// 0-2:keyfrag, 3-8:level, 9:is_topnode  
	int * restrict DUp;				// # nodes till father 
	int * restrict Npart;			// # particles
	
	Float * restrict CoM[3];		// Center of Mass of node
	Float * restrict Dp[3];			// momentum of node
	Float * restrict Rcrit;			// critical radius (Dehnen 2002, eq 10)
#ifdef FMM_SAVE_NODE_POS
	Float * restrict Pos[3];		// center of node
#endif
	union { 						// point to next sibling or leaf
		int * restrict DNext; 		// # nodes to next sibling, if >0
		int * restrict Leaf_Ptr;	// -leaf_idx - 1 ,  if < 0 
	};
	union {
		Float * restrict Mass;		// mass of node is first multipole
		Float * restrict M[NM];		// Multipole of node, eq. 54
	};
	union {
		Float * restrict Force[3];
		Float * restrict Fn[NFN];	// FMM field tensor, eq. 58
	};
} FMM;

size_t sizeof_FMM; // count sizeof of all members

int NNodes, Max_Nodes;

int * restrict Leaf2Part; // convert leaf index to index of first particle
int * restrict Leaf2Node; // convert leaf index to index of containing node

int NLeafs;

omp_lock_t NNodes_Lock, NLeafs_Lock;

double Epsilon[NPARTYPE], // softening
	   Epsilon2[NPARTYPE], 
	   Epsilon3[NPARTYPE];

bool Is_Top_Node(const struct FMM_Node fmm, const int node);
int Level(const struct FMM_Node fmm, const int node);
int Triplet(const struct FMM_Node fmm, const int node);
void Free_Gravity_FMM(struct FMM_Node);
struct FMM_Node Alloc_FMM_Nodes(const int N);
struct FMM_Node Point_FMM_Nodes(const int i);

#else 

static inline void Setup_Gravity_FMM() {};
static inline void Gravity_Acceleration_FMM() {};
static inline void Free_Gravity_FMM(struct FMM_Node) {};
static inline void Gravity_FMM_Build() {};
static inline void Gravity_FMM_P2L() {};
static inline void Gravity_FMM_Update_Kicks() {};
static inline void Gravity_FMM_Update_Topnode_Kicks() {};
static inline void Gravity_FMM_Update_Drift(const double dt) {};

#endif // GRAVITY && GRAVITY_FMM

#endif // FMM_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
