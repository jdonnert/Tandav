#include "globals.h"
#include "proto.h"

#define PARALLEL_THRES 100000


/*
 * Select the kth element out of an array *data with length ndata.
 * Press et al. 1992
 */

static inline void swap(Float * restrict a, Float * restrict b)
{
	Float tmp = *a;
	*a = *b;
	*b = tmp;

	return;
}

Float Select(const int k, const int ndata, Float * restrict data)
{
	int l = 1;
	int ir = ndata;

	for (;;) {

		if (ir <= l+1) { // small partition

			if ((ir == l+1) && (data[ir] < data[l]))
				swap(&data[l], &data[ir]);

			break; // done
		}

		int mid = (l+ir) >> 1;

		swap(&data[mid], &data[l+1]);

		if (data[l] > data[ir])
			swap(&data[l], &data[ir]);

		if (data[l+1] > data[ir])
			swap(&data[l+1], &data[ir]);

		if (data[l] > data[l+1])
			swap(&data[l], &data[l+1]);

		int i = l+1;
		int j = ir;

		Float a = data[l+1];

		for (;;) {

			do
				i++;
			while (data[i] < a);

			do
				j--;
			while (data[j] > a);

			if (j < i)
				break;

			swap(&data[i], &data[j]);
		}

		data[l+1] = data[j];

		data[j] = a;

		if (j >= k)
			ir = j-1;

		if (j <= k)
			l = i;
	}

	return data[k];
}

/*
 * Find an approximation of the median of *data with length ndata in OpenMP.
 */

static int compare_floats(const void * a, const void *b)
{
	const Float *x = (const Float*)a;
	const Float *y = (const Float*)b;

	return (*x > *y) - (*x < *y);
}

static Float *Results = NULL;

Float Median(const int ndata, Float * restrict data)
{
	if (Sim.NThreads == 1 || ndata < PARALLEL_THRES || 1) {
	
		Float median = 0;

		#pragma omp single copyprivate(median)	
		median = Select(ndata >> 1, ndata, data);
		
		return median;

	}
	
	int nPart = ndata * sizeof(*data) / Task.Buffer_Size + 1;
	nPart = MAX(nPart, 32);

	size_t size = ndata / nPart;
	size_t nBytes = size * sizeof(*data);
	int mid = size >> 1;

	#pragma omp single
	Results = Malloc(nPart * sizeof(*Results), "Results");

	Float *buf = Get_Thread_Safe_Buffer(nBytes);

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < nPart; i++) {

		size_t idx = i * size;

		if (i == nPart-1) { // last one

			size = ndata  - idx;
			nBytes = size * sizeof(*data);
			mid = size >> 1;
		}

		memcpy(&buf[0], &data[idx], nBytes);

		Results[i] = Select(mid, size, buf);

	}


	Float median = 0;

	#pragma omp single copyprivate(median)
	{
		median = Select(nPart >> 1, nPart, Results);

		Free(Results);
	}

	return median;
}

void test_median()
{
	const int N = 100000001;
	Float *arr = Malloc(N*sizeof(Float), "arr");
	Float *arr2 = Malloc(N*sizeof(Float), "arr2");
	Float *arr3 = Malloc(N*sizeof(Float), "arr3");

	for (int i = 0; i < N; i++) {

		arr[i] = erand48(Task.Seed);
		arr2[i] = arr[i];
		arr3[i] = arr[i];
	}

	int kth = N/2;
	Float el = Select(kth, N, arr);

	printf("kth = %d el = %g \n", kth, el);

	Float el2 = Median(N, arr3);

	printf("el2 = %g \n", el2);

	Qsort(Sim.NThreads, arr2, N, sizeof(*arr2), &compare_floats);

	printf("%d %g \n", kth, arr2[kth-1]);

	exit(0);

	return ;
}
