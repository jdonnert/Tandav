#include "fmm.h"

#ifdef GRAVITY_FMM

struct FMM_Node FMM = { NULL }; // holds the trees of each top node
int NNodes = 0, Max_Nodes = 0;

size_t sizeof_FMM = 0;

int NLeafs = 0;
int * restrict Leaf2Part = { NULL }; // holds first particle of Leaf
int * restrict Leaf2Node = { NULL }; // holds FMM node of Leaf

omp_lock_t NNodes_Lock, NLeafs_Lock;

void Acceleration_Gravity_FMM()
{
	//if (Sig.Tree_Update) {

	Gravity_FMM_Build(); // also does P2M

	//Gravity_FMM_M2L();

	//}

	//Gravity_FMM_L2L();

	//Gravity_FMM_L2P();

	//Gravity_FMM_P2P();

	return ;
}

void Setup_Gravity_FMM()
{
	omp_init_lock(&NNodes_Lock);
	omp_init_lock(&NLeafs_Lock);

	Max_Nodes = 0.3 * Task.Npart_Total;

	Leaf2Part = Malloc(Task.Npart_Total * sizeof(*Leaf2Part), "Leaf2Part");
	Leaf2Node = Malloc(Task.Npart_Total * sizeof(*Leaf2Node), "Leaf2Node");

	for (int i = 0; i < NPARTYPE; i++) { // Plummer eqiv. softening

		Epsilon[i] = -41.0/32.0 * Param.Grav_Softening[i]; // for Dehnen K1
		Epsilon2[i] = Epsilon[i] * Epsilon[i];
		Epsilon3[i] = Epsilon[i] * Epsilon[i] * Epsilon[i];
	}

	sizeof_FMM += sizeof(FMM.DNext[0]) + sizeof(FMM.Bitfield[0]) 
				+ sizeof(FMM.DUp[0]) + sizeof(FMM.Npart[0]) 
				+ sizeof(FMM.Mass[0]) 
				+ 3*sizeof(FMM.CoM[0][0]) + 3*sizeof(FMM.Dp[0][0]) 
#ifdef FMM_SAVE_NODE_POS
				+ 3*sizeof(FMM.Pos[0][0])
#endif
				+ 0;

	rprintf("\nsizeof(FMM) = %zu byte\n\n", sizeof_FMM);

	return ;
}

/*
 * Handle allocation, deallocation of fmm nodes.
 */

void Gravity_FMM_Free(struct FMM_Node f)
{
	Free(f.DNext);
	Free(f.Bitfield);
	Free(f.DUp);
	Free(f.Npart);
	Free(f.Mass);

	for (int i = 0; i < 3; i++) {

		Free(f.CoM[i]);
		Free(f.Dp[i]);

#ifdef FMM_SAVE_NODE_POS
		Free(f.Pos[i]);
#endif
	}

	return ;
}

struct FMM_Node Alloc_FMM_Nodes(const int N)
{
	struct FMM_Node f = { NULL }; // a collection of pointers

	f.DNext = Malloc(N*sizeof(f.DNext), "FMM.DNext");
	f.Bitfield = Malloc(N*sizeof(f.Bitfield), "FMM.Bitfield");
	f.DUp = Malloc(N*sizeof(f.DUp), "FMM.DUp");
	f.Npart = Malloc(N*sizeof(f.Npart), "FMM.Npart");
	f.Mass = Malloc(N*sizeof(f.Mass), "FMM.Mass");

	f.CoM[0] = Malloc(N*sizeof(f.CoM[0]), "FMM.CoM[0]");
	f.CoM[1] = Malloc(N*sizeof(f.CoM[1]), "FMM.CoM[1]");
	f.CoM[2] = Malloc(N*sizeof(f.CoM[2]), "FMM.CoM[2]");
	f.Dp[0] = Malloc(N*sizeof(f.Dp[0]), "FMM.Dp[0]");
	f.Dp[1] = Malloc(N*sizeof(f.Dp[1]), "FMM.Dp[1]");
	f.Dp[2] = Malloc(N*sizeof(f.Dp[2]), "FMM.Dp[2]");

#ifdef FMM_SAVE_NODE_POS
	f.Pos[0] = Malloc(N*sizeof(f.Pos[0]), "FMM.Pos[0]");
	f.Pos[1] = Malloc(N*sizeof(f.Pos[1]), "FMM.Pos[1]");
	f.Pos[2] = Malloc(N*sizeof(f.Pos[2]), "FMM.Pos[2]");
#endif

	return f;
}

bool Is_Top_Node(const struct FMM_Node fmm, const int node)
{
	return (bool) ((fmm.Bitfield[node] & (1 << 9)) >> 9);
}

#endif // GRAVITY_FMM
