#include "fmm.h"

#ifdef GRAVITY_FMM

static void malloc_nodes(const int N);

uint64_t NNodes = 0, Max_Nodes = 0;
extern struct FMM_Node FMM = { NULL };
static struct FMM_Nodes fmm = { NULL };
static omp_lock_t Node_Lock;

/*
 * This builds the FMM tree in parallel. We OpenMP decompose along top-nodes
 * and vectorize the leafs.
 */

void Gravity_FMM_Build()
{


	return ;
}


void Gravity_FMM_Setup()
{
	omp_init_lock(Node_Lock);

	Max_Nodes = 0.3 * Task.Npart_Total;

	realloc_nodes(Max_Nodes);

	for (int i = 0; i < NPARTYPE; i++) { // Plummer eqiv. softening
	
		Epsilon[i] = -41.0/32.0 * Param.Grav_Softening[i]; // for Dehnen K1
		Epsilon2[i] = Epsilon[i] * Epsilon[i];
		Epsilon3[i] = Epsilon[i] * Epsilon[i] * Epsilon[i];
	}

	return ;
}

static void Free_Nodes()
{
	Free(FMM.DNext);
	Free(FMM.Bitfield);
	Free(FMM.DUp);
	Free(FMM.Npart);
	Free(FMM.Pos[0]);
	Free(FMM.Pos[1]);
	Free(FMM.Pos[2]);
	Free(FMM.Mass);
	Free(FMM.CoM[0]);
	Free(FMM.CoM[1]);
	Free(FMM.CoM[2]);
	Free(FMM.Dp[0]);
	Free(FMM.Dp[1]);
	Free(FMM.Dp[2]);

	return ;
}

static void realloc_nodes(const int N)
{
	if (FMM.DNext == NULL) 
		Free_Nodes();

	FMM.DNext = Malloc(N*sizeof(*FMM.DNext), "FMM.DNext");
	FMM.Bitfield = Malloc(N*sizeof(*FMM.Bitfield), "FMM.Bitfield");
	FMM.DUp = Malloc(N*sizeof(*FMM.DUp), "FMM.DUp");
	FMM.Npart = Malloc(N*sizeof(*FMM.Npart), "FMM.Npart");
	FMM.Pos[0] = Malloc(N*sizeof(*FMM.Pos[0]), "FMM.Pos[0]");
	FMM.Pos[1] = Malloc(N*sizeof(*FMM.Pos[1]), "FMM.Pos[1]");
	FMM.Pos[2] = Malloc(N*sizeof(*FMM.Pos[2]), "FMM.Pos[2]");
	FMM.Mass = Malloc(N*sizeof(*FMM.Mass), "FMM.Mass");
	FMM.CoM[0] = Malloc(N*sizeof(*FMM.CoM[0]), "FMM.CoM[0]");
	FMM.CoM[1] = Malloc(N*sizeof(*FMM.CoM[1]), "FMM.CoM[1]");
	FMM.CoM[2] = Malloc(N*sizeof(*FMM.CoM[2]), "FMM.CoM[2]");
	FMM.Dp[0] = Malloc(N*sizeof(*FMM.Dp[0]), "FMM.Dp[0]");
	FMM.Dp[1] = Malloc(N*sizeof(*FMM.Dp[1]), "FMM.Dp[1]");
	FMM.Dp[2] = Malloc(N*sizeof(*FMM.Dp[2]), "FMM.Dp[2]");
	
	return ;
}

#endif // GRAVITY_FMM
