#include "globals.h"
#include "domain.h"

static int compare_leafs_by_key(const void *a, const void *b);
static void sort_vectors_by_peano_key();

/*
 * We find all leafs in local topnodes using the PH keys. This is basically
 * a tree walk on the bit level ! VECTOR_SIZE sets the max number of 
 * particles in a leaf. Hence we find the smallest level, at which a tree
 * node will contain this many particles.
 * This works by moving a bitmask to select a triplet across the particles 
 * and checking for changes in the PH triplets.
 * */

extern struct Vector_Data Vec = { NULL }; 

void Find_Vectors()
{
	Profile("Find Vec");

	int *first = Get_Thread_Safe_Buffer(Task.Npart_Total_Max*sizeof(*first));

	#pragma omp single // for schedule(static,1) nowait
	for (int i = 1; i < NTop_Nodes; i++) {

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

			for (;;) { // loop over lvl
	
				npart = 1;

				const peanoKey itriplet = P.Key[ipart] & mask;

				for (int jpart = ipart+1; jpart < jmax; jpart++) {

					peanoKey jtriplet = P.Key[jpart] & mask;

					if (jtriplet != itriplet)
						break;

					npart++;
				}

				if (npart <= VECTOR_SIZE) // leaf found
					break;

				mask <<= 3;
				
				lvl++;

				Assert(lvl < 30, "%d %d %d",lvl, ipart, i);
			}
			
			ipart += npart;

			lvl = 0; // find next lvl 
			mask = 0x7;

			while ((P.Key[ipart] & mask) == (P.Key[ipart-1] & mask)) {

				lvl++;
				mask <<= 3;
			}

		} // while ipart

	int start = 0;

	printf("%d %d %d %d \n ", i, NVec, nVec, first[0]);
	#pragma omp critical
	{
		start = NVec;
		NVec += nVec;

	} // omp critical
	
	printf("%d %d \n ", i, NVec);
	for (int j = 0; j > nVec; j++) {
	
		int idx = start + j;
		
		Vec.First[idx] = first[j];
	}

	} // for i

	Vec.First[NVec] = Task.Npart_Total; // add for next loop

	#pragma omp for
	for (int i = 0; i < NVec; i++) {

		Vec.Last[i] = Vec.First[i+1] - 1;
		Vec.Key[i] = P.Key[Vec.First[i]];
	}
	
	sort_vectors_by_peano_key();

	int i = 0;
	printf("%d %d ", Vec.First[i], Vec.Last[i]);
	Print_Int_Bits(Vec.Key[i], 128,0);
	exit(0);

	Profile("Find Vec");

	return ;
}

void Setup_Vectors()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*Vec.First);
	
	Vec.First = Malloc(nBytes, "Vec.First");
	Vec.Last = Malloc(nBytes, "Vec.Last");

	nBytes = Task.Npart_Total_Max*sizeof(*Vec.Key);

	Vec.Key = Malloc(nBytes, "Vec.Key");

	return;
}

static int compare_leafs_by_key(const void *a, const void *b)
{
	const struct Vector_Data *x = (const struct Vector_Data *) a;
	const struct Vector_Data *y = (const struct Vector_Data *) b;
	
	return (int) (x->Key < y->Key) - (x->Key > y->Key);
}

static size_t *Idx = NULL;

static void sort_vectors_by_peano_key()
{
	size_t nBytes = NVec * sizeof(*Idx);

	#pragma omp single
	Idx = Malloc(nBytes, "Leaf Sort Idx");

	Qsort_Index(NThreads, Idx, Vec.Key, NVec, sizeof(*Vec.Key), 
				&compare_leafs_by_key);	

	size_t *idx = Get_Thread_Safe_Buffer(nBytes);

	#pragma omp sections
	{

	#pragma omp section // #0
	{
	
	memcpy(idx, Idx, nBytes);

	Reorder_Array_Char(sizeof(Vec.Key), NVec, Vec.Key, idx);

	} // omp section
	
	#pragma omp section // #1
	{
	
	memcpy(idx, Idx, nBytes);

	Reorder_Array_4(NVec, Vec.First, idx);
	
	} // omp section
	
	#pragma omp section // #2
	{
	
	memcpy(idx, Idx, nBytes);

	Reorder_Array_4(NVec, Vec.Last, idx);

	} // omp section

	} // omp sections
	
	Free(Idx);

	return ;
}
