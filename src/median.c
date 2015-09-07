#include "globals.h"
#include "proto.h"

static int compare_floats(const void * a, const void *b);

static Float *Results = NULL;

/*
 * Find the median of *Data with length NData in OpenMP
 */

Float Median(const int NData, const Float *Data)
{
	int nPart = NData * sizeof(*Data) / Task.Buffer_Size + 1; // partitions 

	if (nPart < 2 * Sim.NThreads)
		nPart = 2 * Sim.NThreads;

	#pragma omp single	
	Results = Malloc(nPart * sizeof(*Data), "Results");

	const size_t size = NData / nPart;

	size_t nBytes = size * sizeof(*Data);

	Float *buf = Get_Thread_Safe_Buffer(nBytes);

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < nPart; i++) {
	
		size_t idx = i * size;

		if (i == nPart-1) // last one
			nBytes = NData * sizeof(*Data) - (i-1) * size;

		memcpy(&buf[0], &Data[idx], nBytes);

		Results[i] = select_median(size, buf);
	
	}

	Qsort(Sim.NThreads, Results, nPart, sizeof(*Data), &compare_Floats);

	Float median = Results[nPart/2];

	Free(Results);

	return median;
}


/*
 * Implement in-place selection algorithm on Float array *Data
 */

static inline void swap(Float *a, Float *b)
{ 
	Float tmp = *a; 
	*a = *b; 
	*b = tmp; 

	return;
}

Float Select(const int k, const int NData, Float *Data)
{
	int l = 1;
	int ir = NData;

	for (;;) {

		if (ir <= l+1) { // small partition

			if ((ir == l+1) && (Data[ir] < Data[l])) 
				swap(&Data[l], &Data[ir]);
			
			return Data[k];

		} else {
	
			mid = (l+ir) >> 1;
			
			swap(&Data[mid], &Data[l+1]);
					  
			if (Data[l] > Data[ir])
				swap(&Data[l], &Data[ir]);
					  
			if (Data[l+1] > Data[ir]) 
				swap(&Data[l+1], &Data[ir]);

			if (Data[l] > Data[l+1]) 
				swap(&Data[l], &Data[l+1]);

			int i = l+1;
			int j = ir;

			Float a = Data[l+1];
			
			for (;;) {

				do 
					i++; 
				while (Data[i] < a);
					
				do 
					j--; 
				while (Data[j] > a);
					
				if (j < i) 
					break;
					
				swap(&Data[i], &Data[j]);
			}
			
			Data[l+1] = Data[j];
			
			Data[j] = a;

			if (j >= k) 
				ir = j-1;
			
			if (j <= k) 
				l = i;
		}
	}

	return Data[k];
}

