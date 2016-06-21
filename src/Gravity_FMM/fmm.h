#ifndef FMM_H
#define FMM_H

#include "../includes.h"

/*
 * Gravity using the Fast Multipole Method with dual tree traversal.
 */

#if defined(GRAVITY) && defined(GRAVITY_FMM)

void Setup_Gravity_FMM();
void Free_Gravity_FMM();
void Gravity_FMM_Build();
void Gravity_FMM_P2L();
void Gravity_FMM_Update_Kicks();
void Gravity_FMM_Update_Topnode_Kicks();
void Gravity_FMM_Update_Drift(const double dt);

extern struct FMM_Node {
	int * restrict DNext;
	uint32_t * restrict Bitfield;
	int * restrict DUp;
	int * restrict Npart;
	Float * restrict Pos[3];
	Float * restrict Mass;
	Float * restrict CoM[3];
	Float * restrict Dp[3];
} FMM;

uint32_t NNodes;

double Epsilon[NPARTYPE], // softening
	   Epsilon2[NPARTYPE], 
	   Epsilon3[NPARTYPE];

int NLeafs;
int * restrict leaf2part;
int * restrict leaf2node;

#else 

inline void Setup_Gravity_FMM() {};
inline void Free_Gravity_FMM() {};
inline void Gravity_FMM_Build() {};
inline void Gravity_FMM_P2L() {};
inline void Gravity_FMM_Update_Kicks() {};
inline void Gravity_FMM_Update_Topnode_Kicks() {};
inline void Gravity_FMM_Update_Drift(const double dt) {};

#endif // GRAVITY && GRAVITY_FMM

#endif // FMM_H
