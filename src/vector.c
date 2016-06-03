#include "vector.h"
#include "Gravity/tree.h"

static int compare_vectors_by_key(const void *a, const void *b);

/*
 * We find all leafs in local topnodes using the PH keys. This is basically
 * a tree walk on the bit level. VECTOR_SIZE sets the max number of 
 * particles in a leaf. Hence we find the smallest level, at which a tree
 * node will contain this many particles.
 * We are using a bitmask to select a triplet across the particles and check
 * for changes in the PH triplets.
 * */

int * restrict Vec = { NULL }; 
int NVec = 0;

static int sum = 0;

void Find_Leaf_Vectors()
{
	Profile("Find Vec");

	int *first = Get_Thread_Safe_Buffer(Task.Npart_Total_Max*sizeof(*first));

	#pragma omp single // for schedule(static,1) nowait
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
	
printf("start=%d NVec=%d nVec=%d \n ", dest, NVec, nVec);

		memcpy(&Vec[dest], first, nVec*sizeof(*first));

	} // for i

	Vec[NVec] = Task.Npart_Total; // add for loops

	Qsort(NThreads, Vec, NVec, sizeof(*Vec), &compare_vectors_by_key);	

	#pragma omp for reduction(+:sum)
	for (int i = 0; i < NVec; i++)
		sum += Vec[i] - Vec[i+1];

	rprintf("Found %d vectors; Average length %d, VECTOR_SIZE %d", 
			NVec, (int)(((double) sum)/NVec), VECTOR_SIZE);

	Profile("Find Vec");
exit(0);
	return ;
}

void Setup_Leaf_Vectors()
{
	size_t nBytes = Task.Npart_Total_Max*sizeof(*Vec);
	
	Vec = Malloc(nBytes, "Vec");

	return;
}

static int compare_vectors_by_key(const void *a, const void *b)
{
	const int i = *((const int *) a);
	const int j = *((const int *) b);
	
	return (int) (P.Key[i] < P.Key[j]) - (P.Key[i] > P.Key[j]);
}

static void test_vectors()
{
	int node = 0;
	int ivec = 0;

	while (node < NNodes) {

		while (Tree[node].DNext > 0)
			node++;
	
		rprintf("Leaf: %d @ %d vec=%d npart=%d ",
			node, Level(node), -Tree[node].DNext-1, Tree[node].Npart );
		rprintf("Vec: vec=%d npart =%d\n", Vec[ivec], 
				Vec[ivec+1] - Vec[ivec]);

		if (ivec == NVec-2) {

		rprintf("Leaf: %d lvl %d vec %d npart %d ",
			node, Level(node), -Tree[node].DNext-1, Tree[node].Npart );
		rprintf("Vec: vec=%d npart =%d\n", Vec[ivec], 
				Vec[ivec+1] - Vec[ivec]);

		int ipart = 31;

		rprintf("%d ", ipart-1); Print_Int_Bits(P.Key[ipart-1], 128,0);
		rprintf("%d ", ipart);  Print_Int_Bits(P.Key[ipart], 128,0);
		rprintf("%d ", ipart+1); Print_Int_Bits(P.Key[ipart+1], 128,0);
		rprintf("%d ", ipart+2); Print_Int_Bits(P.Key[ipart+2], 128,0);
		rprintf("%d ", ipart+3); Print_Int_Bits(P.Key[ipart+3], 128,0);
		rprintf("%d ", ipart+4); Print_Int_Bits(P.Key[ipart+4], 128,0);
		exit(0);
		}

		node++;
		ivec++;
	}
	return ;
}
