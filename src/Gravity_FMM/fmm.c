#include "fmm.h"

#ifdef GRAVITY_FMM

int NLeafs = 0;
int * restrict Leaf2Part = { NULL }; // get first particle of Leaf
int * restrict Leaf2Node = { NULL }; // get FMM node of leaf

void Gravity_Acceleration() 
{
	//if (Sig.Tree_Update) {
	
	find_leaf_vectors();

	Gravity_FMM_Build(); 
	
	Gravity_FMM_P2M();

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

	Leaf_Nodes = Malloc(Task.Npart_Total * sizeof(*Leaf_Node), "Leaf Nodes");

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


/*
 * We find all leafs i local topnodes using the PH keys. This is basically
 * a tree walk on the bit level. VECTOR_SIZE sets the max number of 
 * particles in a leaf. Hence we find the smallest level, at which a tree
 * node will contain this many particles.
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

			if (ipart == last_part) // last particle -> single leaf
				break;

			int jmax = MIN(last_part, ipart + VECTOR_SIZE) + 1;

			int npart = 0; 

			peanoKey mask = ((peanoKey) 0x7) << (3*lvl);

			for (;;) { // loop over level
	
				npart = 1; // account for ipart

				const peanoKey triplet_i = P.Key[ipart] & mask;

				for (int jpart = ipart+1; jpart < jmax; jpart++) {

					peanoKey triplet_j = P.Key[jpart] & mask;

					if (triplet_j != triplet_i)
						break;	

					npart++;
				}

				if (npart <= VECTOR_SIZE) // leaf found
					break;

				mask <<= 3;
				
				lvl++;

				//Assert(lvl < 30, "%d %d %d", lvl, ipart, i);
			}
			
			ipart += npart;

			lvl = D[i].TNode.Level + 1; // find next lvl 
			mask = ((peanoKey) 0x7) << (3*lvl);

			while ((P.Key[ipart] & mask) == (P.Key[ipart-1] & mask)) {

				lvl++;
				mask <<= 3;
			}

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
