#include "fmm.h"

#ifdef GRAVITY_FMM

int NLeafs = 0;
int * restrict Leaf2Part = { NULL }; // holds first particle of Leaf
int * restrict Leaf2Node = { NULL }; // holds FMM node of Leaf

void Gravity_Acceleration() 
{
	//if (Sig.Tree_Update) {
	
	find_leaf_vectors();

	Gravity_FMM_Build_P2M(); 
	
	//Gravity_FMM_M2L();
	
	//} 
	
	//Gravity_FMM_L2L();
	
	//Gravity_FMM_L2P();
	
	//Gravity_FMM_P2P();

	return ;
}

void Setup_Gravity_FMM()
{
	omp_init_lock(Node_Lock);

	Max_Nodes = 0.3 * Task.Npart_Total;

	realloc_nodes(Max_Nodes, FMM);

	Leaf2Part = Malloc(Task.Npart_Total * sizeof(*Leaf_Node), "Leaf Nodes");
	Leaf2Node = Malloc(Task.Npart_Total * sizeof(*Leaf_Node), "Leaf Nodes");

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

void Alloc_Gravity_FMM(const int N, struct FMM_Node *f)
{
	if (f.DNext != NULL) 
		Free_Gravity_FMM(f);

	f->DNext = Malloc(N*sizeof(*f.DNext), "FMM.DNext");
	f->Bitfield = Malloc(N*sizeof(*f.Bitfield), "FMM.Bitfield");
	f->DUp = Malloc(N*sizeof(*f.DUp), "FMM.DUp");
	f->Npart = Malloc(N*sizeof(*f.Npart), "FMM.Npart");
	f->Pos[0] = Malloc(N*sizeof(*f.Pos[0]), "FMM.Pos[0]");
	f->Pos[1] = Malloc(N*sizeof(*f.Pos[1]), "FMM.Pos[1]");
	f->Pos[2] = Malloc(N*sizeof(*f.Pos[2]), "FMM.Pos[2]");
	f->Mass = Malloc(N*sizeof(*f.Mass), "FMM.Mass");
	f->CoM[0] = Malloc(N*sizeof(*f.CoM[0]), "FMM.CoM[0]");
	f->CoM[1] = Malloc(N*sizeof(*f.CoM[1]), "FMM.CoM[1]");
	f->CoM[2] = Malloc(N*sizeof(*f.CoM[2]), "FMM.CoM[2]");
	f->Dp[0] = Malloc(N*sizeof(*f.Dp[0]), "FMM.Dp[0]");
	f->Dp[1] = Malloc(N*sizeof(*f.Dp[1]), "FMM.Dp[1]");
	f->Dp[2] = Malloc(N*sizeof(*f.Dp[2]), "FMM.Dp[2]");

	return ;
}

/*
 * We find all leafs in the local topnodes using the PH keys. This is  a 
 * tree walk on the bit level. VECTOR_SIZE sets the max number of 
 * particles in a leaf. Hence we find the smallest level, at which a tree
 * node will contain up to this many particles.
 * We are using a bitmask to select a triplet across the particles and check
 * for changes in the PH triplets. Obviously particles have to be PH ordered.
 * */

static int compare_vectors(const void *a, const void *b);

static double sum = 0;

void Find_Leaf_Vectors()
{
	Profile("Find Vec");

	int *leafs = Get_Thread_Safe_Buffer(Task.Npart_Total_Max*sizeof(*leafs));

	#pragma omp single
	NLeafs = sum = 0;

	#pragma omp for schedule(static,1) nowait
	for (int i = 0; i < NTop_Nodes; i++) {

		int first_part = D[i].TNode.First_Part;
		int last_part = first_part + D[i].TNode.Npart - 1;

		int n = 0;
		int ipart = first_part;

		int lvl = D[i].TNode.Level + 1;
		
		while (ipart <= last_part) { // all particles in the top node
			
			leafs[n] = ipart;

			n++;

						
			ipart += npart;

			

		} // while ipart

		int dest = 0;

		#pragma omp critical
		{

		dest = NLeafs;
		NLeafs += n;

		} // omp critical
	
		memcpy(&Leaf2Part[dest], leafs, n*sizeof(*leafs));

	} // for i

	Leaf2Part[NVec] = Task.Npart_Total; // terminate for loops

	Qsort(Leaf2Part, NLeafs, sizeof(*Leaf2Part), &compare_leafs);	

	#pragma omp for reduction(+:sum)
	for (int i = 0; i < NVec; i++) 
		sum += Leaf2Part[i+1] - Leaf2Part[i];

	rprintf("Found %d leaf vectors; Average length=%4f, VECTOR_SIZE=%d\n\n", 
			NLeafs, sum/NLeafs, VECTOR_SIZE);

	Profile("Find FMM Leafs");

	return ;
}

void Setup_Leaf_Vectors()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*Vec);
	
	Vec = Malloc(nBytes, "Leaf Vectors");

	return;
}

static int compare_leafs(const void *a, const void *b)
{
	const int i = *((const int *) a);
	const int j = *((const int *) b);
	
	return (int) (Vec[i] < Vec[j]) - (Vec[i] > Vec[j]);
}



#endif // GRAVITY_FMM
