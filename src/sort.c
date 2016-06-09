/* 
 * Simple parallel quicksort, check wikipedia
 * Shamelessly hacked from glibc, thereby GPL2 
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience; 
 * Vol. 23 (11), 1249-1265, 1993. 
 */

#include "sort.h"

#define PARALLEL_THRES_QSORT 15000 // use serial sort below this limit
#define PARALLEL_THRES_HEAPSORT 15000
#define INSERT_THRES 8 // insertion sort threshold



#define PARALLEL_THRESHOLD (1 << 10) // minimum size to use OpenMP
#define MIN_LIB_THRESHOLD (1 << 16) // partition size to switch to std qsort
#define N_PARTITIONS_PER_CPU 16 // number of sub-partition per thread by qsort
#define N_MEDIAN 32 // get pivot element from a median of this

static size_t Lib_Threshold = 0;

double *x, *y;
size_t *p, *q;

static inline void swapAll(void * restrict a, void * restrict b, size_t nBytes)
{	
	char * restrict x = (char *) a;
	char * restrict y = (char *) b;

	char tmp[nBytes];

	memcpy(tmp, x, nBytes);
	memcpy(x, y, nBytes);
	memcpy(y, tmp, nBytes);

	return ;
}

static void swap1(void * restrict a, void * restrict b)
{
	char * restrict x = (char *) a;
	char * restrict y = (char *) b;

	char tmp = *x;

	*x = *y;
	*y = tmp;

	return ;
}

static void swap2(void * restrict a, void * restrict b)
{
	uint16_t * restrict x = (uint16_t *) a;
	uint16_t * restrict y = (uint16_t *) b;

	uint16_t tmp = *x;

	*x = *y;
	*y = tmp;

	return ;
}

static void swap4(void * restrict a, void * restrict b)
{
	uint32_t * restrict x = (uint32_t *) a;
	uint32_t * restrict y = (uint32_t *) b;

	uint32_t tmp = *x;

	*x = *y;
	*y = tmp;

	return ;
}

static void swap8(void * restrict a, void * restrict b)
{
	uint64_t * restrict x = (uint64_t *) a;
	uint64_t * restrict y = (uint64_t *) b;

	uint64_t tmp = *x;

	*x = *y;
	*y = tmp;

	return ;
}

static void swap16(void * restrict a, void * restrict b)
{
	__uint128_t * restrict x = (__uint128_t *) a;
	__uint128_t * restrict y = (__uint128_t *) b;

	__uint128_t tmp = *x;

	*x = *y;
	*y = tmp;

	return ;
}

static void (*swap) ();

static inline void swap_size_t(size_t * restrict a, size_t * restrict b) 
{
	size_t tmp = *a;

	*a = *b;
	*b = tmp;
		
	return ;
}

#define COMPARE_DATA(a,b,size) ((*cmp) ((char *)(data + *a * size), \
	(char *)(data + *b * size))) // compare with non-permutated data

static char * median_of(const int N, char *lo, size_t nData, size_t size,
						int (*cmp)  (const void*, const void *))
{
	char *hi = lo + (nData-1) * size;
	size_t dp = (nData / N) * size;

	char *addr[N];

	for (int i = 0; i < N; i++)
		addr[i] = lo + i*dp;
	addr[N-1] = hi;

	for (int i = 1; i < N; i++) // insertion sort
		for (int j = i; j > 0 && cmp(addr[j-1],addr[j])>0; j--)
			swap(addr[j], addr[j-1]);

	return addr[N>>1];
}

static void omp_qsort(void *Data, size_t nData, size_t size, 
			int (*cmp) (const void*, const void *))
{
	char *lo = Data;
	char *hi = Data + (nData-1)*size;

	char *mid = median_of(N_MEDIAN, lo, nData, size, cmp); 

	char *left  = lo + size;
	char *right = hi - size;
	 
	do { // partition

		while ((*cmp)((void *)left, (void *)mid) < 0)
			left += size;

		while ((*cmp)((void *)mid,(void *)right) < 0)
			right -= size;
		
		if (left < right) {

			(*swap)(left, right);

			if (mid == left)
				mid = right;
			else if (mid == right)
				mid = left;

			left += size;
			right -= size;

		} else if (left == right) {

			left += size;
			right -= size;

			break;
		}

	} while (left <= right);

	size_t nLeft = (right - lo) / size + 1; // kick off new partitions
	size_t nRight = (hi - left) / size + 1;

	if (nLeft > 1) {
	
		if (nLeft < Lib_Threshold) {

			#pragma omp task 
			qsort(lo, nLeft, size, cmp); // qsort is likely pretty good ...

		} else {
		
			#pragma omp task 
			omp_qsort(lo, nLeft, size, cmp);
		
		}
	}

	if (nRight > 1) {
	
		if (nRight < Lib_Threshold) {

			#pragma omp task 
			qsort(left, nRight, size, cmp);

		} else {

			#pragma omp task 
			omp_qsort(left, nRight, size, cmp);
		}
	}

	return ;
}

/*
 * A thread safe OpenMP quicksort. We use a function pointer for the optimised
 * swap routine. This is not as bad as it sounds with optimised compilers. We
 * rely on the build-in qsort whereever we can, because it will be highly
 * optimised. For large arrays we find the real median of a few values using
 * a insertion sort and use it as a pivot. (Bentley & McIlroy 1993)
 */

void Qsort(void *Data, size_t nData, size_t size, 
			int (*cmp) (const void*, const void *))
{
	if ( (nData < PARALLEL_THRESHOLD) || (! omp_in_parallel()) ) {
		
		#pragma omp single
		qsort(Data, nData, size, cmp);

		return ;
	}

	switch (size) { // set swap function pointer
	
		case 1: swap = &swap1;
				break;

		case 2: swap = &swap2;
				break;

		case 4: swap = &swap4;
				break;
		
		case 8: swap = &swap8;
				break;
		
		case 16: swap = &swap16;
				break;
	}

	Lib_Threshold = nData / NThreads / N_PARTITIONS_PER_CPU; // load balancing
	Lib_Threshold = MAX(PARALLEL_THRESHOLD, Lib_Threshold);
	
	#pragma omp single
	omp_qsort(Data, nData, size, cmp); // burn baby !
	
	return ;
}
/* testing */
int test_compare(const void * a, const void *b)
{
	const double *x = (const double*)a;
	const double *y = (const double*)b;

	return (*x > *y) - (*x < *y);
}

void test_sort()
{
	const int Nit = 2;
	int good;

	printf("Testing sort: 1 - 2^31, %d iterations\n"
			"MIN_LIB_THRESHOLD %d \nPARALLEL_THRESHOLD %d\n"
			" N_PARTITIONS_PER_CPU %d \n\n"
			,Nit, MIN_LIB_THRESHOLD, PARALLEL_THRESHOLD, N_PARTITIONS_PER_CPU);

	for (int N = 1ULL << 10; N < (1ULL << 31); N<<=1) {

		x = (double *) malloc( N * sizeof(*x) );
		y = (double *) malloc( N * sizeof(*y) );
	//	p = (size_t *) malloc( N * sizeof(*p) );
	//	q = (size_t *) malloc( N * sizeof(*q) );

		clock_t time = clock(), time2 = clock(), time3 = clock();
		double deltasum0 = 0, deltasum1 = 0;

		rprintf("%3g %10d ", log2(N), N);


	/* in-place sort */

	for (int i = 0; i < Nit; i++) {
		
		#pragma omp parallel
		{
		
		#pragma omp for
		for (int j = 0; j < N; j++) 
    		y[j] = erand48(Task.Seed);

		#pragma omp master
		{

		memcpy(x, y, N*sizeof(*x));
		
		time = clock();

		} // master
			
		#pragma omp barrier
  		
		Qsort(x, N, sizeof(*x), &test_compare);

		#pragma omp master
		{

  		good = 1;

	  	for (int i = 1; i < N; i++)
			if (x[i] < x[i-1])
				good = 0;

		if (good == 0) {
		
			printf("ERROR: Array not sorted :-( \n");

			exit(0);
		}

		time2 = clock();
		deltasum0 += time2-time;

		} // master
		
		#pragma omp barrier

		} // omp parallel

		memcpy(x, y, N* sizeof(*x));

		time = clock();
	
 		qsort(x, N, sizeof(*x), &test_compare);
  	
	 	time2 = clock();

		deltasum1 += time2-time;
	}
		
	deltasum0 /= Nit; 
	deltasum1 /= Nit;

  	printf("%6e %6e %4g \n",
		deltasum0/CLOCKS_PER_SEC/NThreads, 
		deltasum1/CLOCKS_PER_SEC,deltasum1/deltasum0*NThreads );

	free(x); free(p); free(q); free(y);

	} // for N

	exit(0);

	return ;
}

#define STACK_SIZE (CHAR_BIT * sizeof(size_t))
/* 
 * This is an OpenMP parallel external sort.
 * We use the same algorithm as above to divide and conquer. In the first stage
 * we sort 1 partition per thread. The work is inserted in the shared stack
 * with a delta that shrinks from N partitions / 2 to 1. Some partitions can be
 * empty, when a previous partition is too small to split.
 * The second stage is a threaded qsort on the subpartitions down to 
 * INSERT_THRES followed by insertion sort on the nearly ordered subpartition. 
 * Here we are swapping ONLY the permutation array *perm and  are comparing 
 * ONLY to the data array indexed by *perm. I.e. data[] remains unchanged
 */

static struct SharedStackDataSizeT {
    size_t *lo;
    size_t *hi;
} shared_stack_sizet[STACK_SIZE] = {{0,0}};

void Qsort_Index(const int nThreads, size_t *perm, void * const data,
		const int nData, const size_t datasize,
		int (*cmp) (const void *, const void *))
{
	Assert(perm != NULL, "*perm is a NULL pointer, no space to sort");
	Assert(data != NULL, "*data is a NULL pointer, no space to sort");

	if (nData < PARALLEL_THRES_HEAPSORT || nThreads == 1 || true) {

		#pragma omp single
		gsl_heapsort_index(perm, data, nData, datasize, cmp);

		#pragma omp flush

		return ;
	}

	const int desNumPar = MIN(nThreads, floor(nData/INSERT_THRES));

	#pragma omp for
	for (size_t i = 0; i < nData; i++ )
		perm[i] = i;

	#pragma omp flush

	#pragma omp single
	{

	shared_stack_sizet[0].lo = perm;
	shared_stack_sizet[0].hi = &perm[nData - 1];

	} // omp single

	int delta = desNumPar << 1; // twice the distance of entries in shared stack

	/* First stage: subpartitions, roughly NThreads */

	while (delta > 2) {

		delta >>= 1; // we compensated for this before

		#pragma omp flush

		#pragma omp for schedule(static,1)
		for (int i = 0; i < desNumPar; i += delta) {

			size_t *lo = shared_stack_sizet[i].lo;
			size_t *hi = shared_stack_sizet[i].hi;

			if (hi - lo < 2)
				continue;

			size_t *mid = lo + ((hi - lo) >> 1); // pivot & presort

			if (COMPARE_DATA(mid,lo,datasize) < 0) 
				swap_size_t (mid, lo);

			if (COMPARE_DATA(hi,mid,datasize) < 0)
				swap_size_t (hi, mid);
			else
				goto jump_over;

			if (COMPARE_DATA(mid,lo,datasize) < 0) 
				swap_size_t (mid, lo);

			jump_over:;

			size_t *left  = lo + 1;
	 		size_t *right = hi - 1;

		 	do { // collapse the walls 

				while (COMPARE_DATA(left,mid,datasize) < 0) 
					left++;

				while (COMPARE_DATA(mid,right,datasize) < 0)
					right--;

				if (left < right) { 

					swap_size_t (left, right);

					if (mid == left)
			   			mid = right;
					else if (mid == right)
		   				mid = left;
 
					left++;
					right--;

				} else if (left == right) {

					left++; 

					break;
				}

			} while (left <= right);

			shared_stack_sizet[i].lo = lo; // Push to replace old position
			shared_stack_sizet[i].hi = right;

			shared_stack_sizet[i + (delta>>1)].lo = left; //Push for next
			shared_stack_sizet[i + (delta>>1)].hi = hi;
		} // j

    } // while

	/* Second stage, every thread: Qsort on subpartition *beg to *end */

	#pragma omp barrier

	size_t *beg = shared_stack_sizet[Task.Thread_ID].lo;
	size_t *end = shared_stack_sizet[Task.Thread_ID].hi;

	size_t partition_size = end - beg + 1;

	if (partition_size < INSERT_THRES) 
		goto insertion_sort;

	struct SharedStackDataSizeT stack[STACK_SIZE] = { {beg,end} }; 

	int next = 0; 
			
	for (;;)  { 

		size_t *lo = stack[next].lo;  // pop from private stack
		size_t *hi = stack[next--].hi;

  		size_t *mid = lo + ((hi - lo) >> 1); 
		
		if (COMPARE_DATA(mid,lo,datasize) < 0) 
			swap_size_t (mid, lo);

		if (COMPARE_DATA(hi,mid,datasize) < 0)  
			swap_size_t (mid, hi);
		else
	    	goto hop_over;
  	
		if (COMPARE_DATA(mid,lo,datasize) < 0) 
			swap_size_t (mid, lo);
	
		hop_over:;

	  	size_t *left  = lo + 1;
  		size_t *right = hi - 1;

	  	do { // collapse the walls
		
			while (COMPARE_DATA(left,mid,datasize) < 0) 
				left++;

			while (COMPARE_DATA(mid,right,datasize) < 0)
				right--;

    		if (left < right) { 

				swap_size_t (left, right);

				if (mid == left)
		   			mid = right;
				else if (mid == right)
	   				mid = left;
	  
				left++;
				right--;
	
			} else if (left == right) {
	  			
				left++;
				right--;
	 
				break;
			}

		} while (left <= right);

		if (right - lo < hi - left) { // push smaller partition first 

			if (right - lo + 1 > INSERT_THRES) { // push

				stack[++next].lo = lo;
				stack[next].hi = right;
			}

			if (hi - left + 1 > INSERT_THRES) {

				stack[++next].lo = left;
				stack[next].hi = hi;
			}

		} else {

			if (hi - left + 1 > INSERT_THRES) {

				stack[++next].lo = left;
				stack[next].hi = hi;
			}

			if (right - lo + 1 > INSERT_THRES) { // push

				stack[++next].lo = lo;
				stack[next].hi = right;
			}
		}
		
		if (next < 0) // stack empty
			break;

	} // for()
   
	/* Third stage, every thread: insertion sort on its partition */

	insertion_sort:;

	size_t *trail = beg;

	const size_t *runMax = MIN(end, beg + INSERT_THRES);

	for (size_t *run = beg+1; run <= runMax; run++) // smallest element first
		if (COMPARE_DATA(run, trail, datasize) < 0)
			trail = run;

	if (trail != beg)
		swap_size_t(trail, beg);

	size_t *run = beg + 1;

	while (++run <= end) { // insertion sort left to right

		trail = run - 1;

		while (COMPARE_DATA(run, trail, datasize) < 0) // find correct place
			trail--;

		trail++;

		if (trail != run)  { // now move block one forward

			size_t save = *run;

			size_t *hi = run;
			size_t *lo = hi - 1;

			while (lo >= trail) 
				*hi-- = *lo--;

			*hi = save; // sort in *run element
		}
	} // while

	#pragma omp barrier

	return;
}



