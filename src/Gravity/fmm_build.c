#include "fmm.h"

#ifdef GRAVITY_FMM

static void realloc_nodes(const int N);
void Free_Nodes();

uint64_t NNodes = 0, Max_Nodes = 0;
extern struct FMM_Node FMM = { NULL };
static struct FMM_Node fmm = { NULL };
#pragma omp threadprivate(fmm)
static omp_lock_t Node_Lock;

/*
 * This builds the FMM tree in parallel. We OpenMP decompose along top-nodes
 * and vectorize the leafs.
 */

void Gravity_FMM_Build()
{
	Profile("Grav FMM Build");

	#pragma omp single
	NNodes = 0;


	Profile("Grav FMM Build");

	return ;
}


void Gravity_FMM_Setup()
{
	omp_init_lock(Node_Lock);

	Max_Nodes = 0.3 * Task.Npart_Total;

	realloc_nodes(Max_Nodes, FMM);

	for (int i = 0; i < NPARTYPE; i++) { // Plummer eqiv. softening
	
		Epsilon[i] = -41.0/32.0 * Param.Grav_Softening[i]; // for Dehnen K1
		Epsilon2[i] = Epsilon[i] * Epsilon[i];
		Epsilon3[i] = Epsilon[i] * Epsilon[i] * Epsilon[i];
	}

	return ;
}

void Free_Gravity_FMM(struct FMM_Node *f)
{
	Free(f->DNext);
	Free(f->Bitfield);
	Free(f->DUp);
	Free(f->Npart);
	Free(f->Pos[0]);
	Free(f->Pos[1]);
	Free(f->Pos[2]);
	Free(f->Mass);
	Free(f->CoM[0]);
	Free(f->CoM[1]);
	Free(f->CoM[2]);
	Free(f->Dp[0]);
	Free(f->Dp[1]);
	Free(f->Dp[2]);

	return ;
}

static void realloc_nodes(const int N, struct FMM_Node *f)
{
	if (f.DNext != NULL) 
		Free_Nodes(f);

	f->DNext = Malloc(N*sizeof(*FMM.DNext), "FMM.DNext");
	f->Bitfield = Malloc(N*sizeof(*FMM.Bitfield), "FMM.Bitfield");
	f->DUp = Malloc(N*sizeof(*FMM.DUp), "FMM.DUp");
	f->Npart = Malloc(N*sizeof(*FMM.Npart), "FMM.Npart");
	f->Pos[0] = Malloc(N*sizeof(*FMM.Pos[0]), "FMM.Pos[0]");
	f->Pos[1] = Malloc(N*sizeof(*FMM.Pos[1]), "FMM.Pos[1]");
	f->Pos[2] = Malloc(N*sizeof(*FMM.Pos[2]), "FMM.Pos[2]");
	f->Mass = Malloc(N*sizeof(*FMM.Mass), "FMM.Mass");
	f->CoM[0] = Malloc(N*sizeof(*FMM.CoM[0]), "FMM.CoM[0]");
	f->CoM[1] = Malloc(N*sizeof(*FMM.CoM[1]), "FMM.CoM[1]");
	f->CoM[2] = Malloc(N*sizeof(*FMM.CoM[2]), "FMM.CoM[2]");
	f->Dp[0] = Malloc(N*sizeof(*FMM.Dp[0]), "FMM.Dp[0]");
	f->Dp[1] = Malloc(N*sizeof(*FMM.Dp[1]), "FMM.Dp[1]");
	f->Dp[2] = Malloc(N*sizeof(*FMM.Dp[2]), "FMM.Dp[2]");
	
	return ;
}

static void prepare_fmm()
{
	if (NNodes != 0) { // tree build aborted, increase mem

		Max_Nodes = TREE_ENLARGEMENT_FACTOR * Max_Nodes;
		
		printf("(%d:%d) Increased tree memory to %6.1f MB, "
			"max %10d nodes, ratio %4g \n", Task.Rank, Task.Thread_ID, 
			Max_Nodes * sizeof(*Tree)/1024.0/1024.0, Max_Nodes, 
			(double) Max_Nodes/Task.Npart_Total); 
	}

	realloc_nodes(Max_Nodes, FMM);
	Tree = Realloc(Tree, Max_Nodes * sizeof(*Tree), "Tree");
		
	memset(Tree, 0, Max_Nodes * sizeof(*Tree));

	NNodes = 0;

	return ;
}

#endif // GRAVITY_FMM
