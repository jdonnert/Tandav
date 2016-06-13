#include "vector.h"
#include "Gravity/tree.h"

static int compare_vectors(const void *a, const void *b);

/*
 * We find all leafs in local topnodes using the PH keys. This is basically
 * a tree walk on the bit level. VECTOR_SIZE sets the max number of 
 * particles in a leaf. Hence we find the smallest level, at which a tree
 * node will contain this many particles.
 * We are using a bitmask to select a triplet across the particles and check
 * for changes in the PH triplets. Obviously particles have to be PH ordered.
 * */

int * restrict Vec = { NULL }; 
int NVec = 0;

static double sum = 0;

void Find_Leaf_Vectors()
{
	Profile("Find Vec");

	int *first = Get_Thread_Safe_Buffer(Task.Npart_Total_Max*sizeof(*first));

	#pragma omp single
	{
	NVec = 0;
	} // omp single

	#pragma omp for schedule(static,1) nowait
	for (int i = 0; i < NTop_Nodes; i++) {

		int first_part = D[i].TNode.First_Part;
		int last_part = first_part + D[i].TNode.Npart - 1;

		int nVec = 0;
		int ipart = first_part;

		int lvl = D[i].TNode.Level + 1;
		
		while (ipart <= last_part) { // all particles in the top node

			first[nVec] = ipart;

			nVec++;

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

		dest = NVec;
		NVec += nVec;

		} // omp critical
	
		memcpy(&Vec[dest], first, nVec*sizeof(*first));

	} // for i

	Vec[NVec] = Task.Npart_Total; // terminate for loops

	Qsort(Vec, NVec, sizeof(*Vec), &compare_vectors);	

	sum = 0;

	#pragma omp for reduction(+:sum)
	for (int i = 0; i < NVec; i++) 
		sum += Vec[i+1] - Vec[i];

	rprintf("Found %d leaf vectors; Average length=%g, VECTOR_SIZE=%d\n\n", 
			NVec, sum/NVec, VECTOR_SIZE);

	Profile("Find Vec");

	return ;
}

void Setup_Leaf_Vectors()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*Vec);
	
	Vec = Malloc(nBytes, "Leaf Vectors");

	return;
}

static int compare_vectors(const void *a, const void *b)
{
	const int i = *((const int *) a);
	const int j = *((const int *) b);
	
	return (int) (Vec[i] < Vec[j]) - (Vec[i] > Vec[j]);
}

