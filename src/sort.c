/* Simple parallel sorting based on subdiving the list */

#include "globals.h"

#define THREADEDLIMIT 4096 // Num of elements below which we use std qsort

void Omp_Sort(char *buf, const size_t num, const size_t size, 
		int (*compare)  (const void *, const void *))
{
	const int nTID = Sim.NThreads;
	const int tID = Task.ThreadID;

	#pragma omp parallel
	if (num > THREADEDLIMIT) {
		
		size_t stride = num*size / nTID + 1;
		
		size_t istart = tID * stride;
		
		for (int i = 1; i < nTID/2; i++) {
		
			if (istart < num)
				qsort(&buf[istart*size], stride, size, compare);

			stride *= 2;
			istart *= 2;

			#pragma omp barrier
		}

	} 

	qsort(buf, num, size, compare);

	return ;
}

#undef THREADEDLIMIT 
