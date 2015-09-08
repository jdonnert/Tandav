#include "globals.h"
#include "proto.h"

#define PARALLEL_THRES 100000

static int compare_floats(const void * a, const void *b);

static Float *Results = NULL;

/*
 * Find an approximation of the median of *data with length ndata in OpenMP.
 */

Float Median(const int ndata, Float *data)
{
	if (Sim.NThreads == 1 || ndata < PARALLEL_THRES)
		return Select(ndata >> 1, ndata, data);

	int nPart = ndata * sizeof(*data) / Task.Buffer_Size + 1; // partitions 

	if (nPart < Sim.NThreads)
		nPart = Sim.NThreads;

	size_t size = ndata / nPart;
	size_t nBytes = size * sizeof(*data);
	int mid = size >> 1;

	#pragma omp single
	Results = Malloc(nPart * sizeof(*data), "Results");

	Float *buf = Get_Thread_Safe_Buffer(nBytes);

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < nPart; i++) {

		size_t idx = i * size;

		if (i == nPart-1) { // last one

			size = (ndata  - idx);
			nBytes = size * sizeof(*data);
			mid = size >> 1;
		}

		memcpy(&buf[0], &data[idx], nBytes);

		Results[i] = Select(mid, size, buf);

	}

	Float median = Select(nPart >> 1, nPart, Results);

	#pragma omp single
	Free(Results);

	return median;
}

static int compare_floats(const void * a, const void *b)
{
	const Float *x = (const Float*)a;
	const Float *y = (const Float*)b;

	return (*x > *y) - (*x < *y);
}

/*
 * Implement in-place selection algorithm on Float array *data ,
 * Press et al. 1992
 */

static inline void swap(Float *a, Float *b)
{
	Float tmp = *a;
	*a = *b;
	*b = tmp;

	return;
}

Float Select(const int k, const int ndata, Float *data)
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
}
