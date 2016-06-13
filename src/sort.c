#include "sort.h"

#define PARALLEL_THRESHOLD (1<<10)	// minimum size to use OpenMP
#define INSERTION_THRESHOLD (1<<5)	// switch to insertion sort here
#define MEDIAN_THRESHOLD (1<<6)		// switch to median of 9 here

#define N_PARTITIONS_PER_CPU 32   	// # sub-partition per thread before qsort

#define CMP_DATA(a,b,size) ((*cmp) ((char *)data + *(a) * size, \
									(char *)data + *(b) * size)) 	

static void set_swap_function(const int size);
static void swap1(void * restrict a, void * restrict b);
static void swap2(void * restrict a, void * restrict b);
static void swap4(void * restrict a, void * restrict b);
static void swap8(void * restrict a, void * restrict b);
static void swap16(void * restrict a, void * restrict b);
static void swap_size_t(size_t * restrict a, size_t * restrict b);

static void omp_qsort(void *data, size_t ndata, size_t size, 
			int (*cmp) (const void*, const void *));
static void omp_qsort_index(size_t *perm, void *data, size_t ndata, size_t size,
			int (*cmp) (const void*, const void *));

static size_t Spawn_Threshold = 0;
static void (*swap) ();

/*
 * A thread safe OpenMP quicksort, in-place variant. 
 * We use a function pointer for the optimised swap routine. This is not as 
 * bad as it sounds with modern compilers. We rely on the build-in qsort 
 * whereever we can, because it is likely very optimised. For large arrays 
 * we use the real median of 9 values as pivot and presort using an insertion 
 * sort. 
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience; 
 * Vol. 23 (11), 1249-1265, 1993.
 */

void Qsort(void *data, size_t ndata, size_t size, 
		   int (*cmp) (const void*, const void *))
{
	Assert(data != NULL, "You gave me a NULL pointer to sort");

	if ( (ndata < PARALLEL_THRESHOLD) || (! omp_in_parallel()) ) {
		
		#pragma omp single
		qsort(data, ndata, size, cmp);

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
	
	Spawn_Threshold = ndata / NThreads / N_PARTITIONS_PER_CPU; // load balancing
	Spawn_Threshold = MAX(PARALLEL_THRESHOLD, Spawn_Threshold);

	#pragma omp single
	omp_qsort(data, ndata, size, cmp); // burn baby !
	
	return ;
}
/*
 * Thread safe OpenMP quicksort, external variant. 
 * We sort the permutation array, but compare the data array at the permuted
 * indices.
 */
void Qsort_Index(size_t *perm, void * data, size_t ndata, size_t size,
				 int (*cmp) (const void *, const void *))
{
	Assert(data != NULL, "You gave me a NULL pointer to sort");
	Assert(perm != NULL, "You gave me a NULL pointer for the permutations");

	#pragma omp for
	for (size_t i = 0; i < ndata; i++)
		perm[i] = i;

	Spawn_Threshold = ndata / NThreads / N_PARTITIONS_PER_CPU * 2; 
	Spawn_Threshold = MAX(PARALLEL_THRESHOLD, Spawn_Threshold);

	#pragma omp single
	omp_qsort_index(perm, data, ndata, size, cmp);
	
	return ;
}

/*
 * Find median and presort
 */

static char * median_of_9(char *lo, size_t ndata, size_t size,
						int (*cmp)  (const void*, const void *))
{
	char *hi = lo + (ndata-1) * size;
	size_t dp = (ndata / 9) * size;

	char *addr[9] = { lo, lo+dp, lo+2*dp, lo+3*dp, lo+4*dp, lo+5*dp, lo+6*dp,
					  lo + 7*dp, hi};

	for (int i = 1; i < 9; i++) // insertion sort
		for (int j = i; j > 0 && cmp(addr[j-1],addr[j]) > 0; j--)
			(*swap)(addr[j], addr[j-1]);

	return addr[5];
}

static void omp_qsort(void *data, size_t ndata, size_t size, 
					  int (*cmp) (const void*, const void *))
{
	char *lo = data;
	char *hi = data + (ndata-1)*size;

	char *mid = median_of_9(lo, ndata, size, cmp); 

	char *left  = lo + size;
	char *right = hi - size;
	 
	do { // partition

		while ((*cmp)((void *)left, (void *)mid) < 0)
			left += size;

		while ((*cmp)((void *)mid, (void *)right) < 0)
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
	
		if (nLeft < Spawn_Threshold) { // stop creating more tasks

			#pragma omp task 
			qsort(lo, nLeft, size, cmp); // qsort is likely pretty good ...

		} else {
		
			#pragma omp task 
			omp_qsort(lo, nLeft, size, cmp);
		
		}
	}

	if (nRight > 1) {
	
		if (nRight < Spawn_Threshold) {

			#pragma omp task 
			qsort(left, nRight, size, cmp);

		} else {

			#pragma omp task 
			omp_qsort(left, nRight, size, cmp);
		}
	}

	return ;
}

static size_t * median_of_9_index(size_t *lo, void *data, size_t ndata, 
		size_t size, int (*cmp)  (const void*, const void *))
{
	size_t *hi = lo + (ndata-1);
	size_t dp = (ndata / 9);

	size_t *addr[9] = { lo, lo+dp, lo+2*dp, lo+3*dp, lo+4*dp, lo+5*dp, lo+6*dp,
					  lo + 7*dp, hi};

	for (int i = 1; i < 9; i++) // insertion sort on perm
		for (int j = i; j > 0 && CMP_DATA(addr[j-1],addr[j], size) > 0; j--)
			swap_size_t(addr[j], addr[j-1]);

	return addr[5]; // return median of 9
}

static size_t * median_of_3_index(size_t *lo, void *data, size_t ndata, 
		size_t size, int (*cmp)  (const void*, const void *))
{
	
	size_t *hi = lo + ndata - 1;
	size_t *mid = lo + (ndata >> 1);

	if (CMP_DATA(mid, lo, size) < 0) 
    	swap_size_t(mid, lo);
  
	if (CMP_DATA(hi, mid, size) < 0)
		swap_size_t(hi, mid);
    else
    	goto jump_over;
  
    if (CMP_DATA(mid, lo, size) < 0) 
    	swap_size_t(mid, lo);
  
    jump_over:;

	return mid;
}


static void omp_qsort_index(size_t *perm, void *data, size_t ndata, size_t size,
							int (*cmp) (const void*, const void *))
{
	if (ndata < INSERTION_THRESHOLD) {

		size_t *lo = perm;
		size_t *hi = lo + ndata;
	
		for (size_t *r = lo + 1; r < hi; r++) // insertion sort
			for (size_t *tr = r; tr > lo && CMP_DATA(tr-1, tr, size) > 0; tr--)
				swap_size_t(tr, tr-1);

		return ;
	}

	size_t *lo = perm;
	size_t *hi = perm + (ndata-1);

	size_t *mid = NULL;
	
	if (ndata < MEDIAN_THRESHOLD)
		mid = median_of_3_index(perm, data, ndata, size, cmp); 
	else
		mid = median_of_9_index(perm, data, ndata, size, cmp); 
	
	size_t *left  = lo + 1;
	size_t *right = hi - 1;
	
	do { // partition

		while (CMP_DATA(left, mid, size) < 0)
			left++;

		while (CMP_DATA(mid, right, size) < 0)
			right--;
		
		if (left < right) {

			swap_size_t(left, right);

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

	size_t nLeft = right - lo + 1; // kick off new partitions
	size_t nRight = hi - left + 1;

	if (nLeft > 1) { // GSL heapsort is so slow, we use our own qsort_index ... 
	
		if (nLeft > Spawn_Threshold) {

			#pragma omp task 
			omp_qsort_index(lo, data, nLeft, size, cmp);

		} else {
		
			omp_qsort_index(lo, data, nLeft, size, cmp);
		}
	}

	if (nRight > 1) {

		if (nRight > Spawn_Threshold) {

			#pragma omp task
			omp_qsort_index(left, data, nRight, size, cmp);
 		
		} else {
	
			omp_qsort_index(left, data, nRight, size, cmp);
		}
	}

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


static inline void swap_size_t(size_t * restrict a, size_t * restrict b) 
{
	size_t tmp = *a;

	*a = *b;
	*b = tmp;
		
	return ;
}


/* Testing */

int test_compare(const void * a, const void *b)
{
	const double *x = (const double*)a;
	const double *y = (const double*)b;

	return (*x > *y) - (*x < *y);
}


static double *x, *y;
static size_t *p, *q;

void test_sort()
{
	const int Nit = 8;
	size_t Nmax = (size_t) 1 << 33;
	int good;

	x = (double *) malloc( Nmax * sizeof(*x) );
	y = (double *) malloc( Nmax * sizeof(*y) );

	printf("Testing sort: 1 -> 2^%g, %d iterations\n"
			"PARALLEL_THRESHOLD = %d\n"
			"N_PARTITIONS_PER_CPU = %d\n", log2(Nmax),
			Nit, PARALLEL_THRESHOLD, N_PARTITIONS_PER_CPU);
	
	for (size_t N = 1ULL << 9; N < Nmax; N<<=1) {

		clock_t time = clock(), time2 = clock(), time3 = clock();
		double deltasum0 = 0, deltasum1 = 0;

		rprintf("%3g %10zu %4.1f GB | ", log2(N), N, N*sizeof(*x)/p3(1024.0));


	for (int i = 0; i < Nit; i++) {  // in-place sort 
		
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

	} // for N
	
	/* external sort */

	printf("Testing external sort: 1 -> 1<<%g, %d iterations\n"
			"INSERTION_THRESHOLD = %d \n"
			"N_PARTITIONS_PER_CPU = %d\n"
			"MEDIAN_THRESHOLD = %d \n", 
			log2(Nmax),
			Nit, INSERTION_THRESHOLD, N_PARTITIONS_PER_CPU, MEDIAN_THRESHOLD);
	
	p = (size_t *) malloc( Nmax * sizeof(*p) );

	for (size_t N = 1ULL << 9; N < Nmax; N<<=1) {
		
		clock_t time = clock(), time2 = clock(), time3 = clock();
		double deltasum0 = 0, deltasum1 = 0;

		rprintf("%3g %10zu %4.1f GB | ", log2(N), N, N*sizeof(*x)/p3(1024.0));
	
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
  		
		Qsort_Index(p, x, N, sizeof(*x), &test_compare);

		#pragma omp master
		{

  		good = 1;

	  	for (int i = 1; i < N; i++)
			if (x[p[i]] < x[p[i-1]])
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
	
	 	gsl_heapsort_index(p, x, N, sizeof(*x), test_compare);
  	
	 	time2 = clock();

		deltasum1 += time2-time;

		deltasum0 /= Nit; 
		deltasum1 /= Nit;

  		printf("%6e %6e %4g \n",
			deltasum0/CLOCKS_PER_SEC/NThreads, 
			deltasum1/CLOCKS_PER_SEC,deltasum1/deltasum0*NThreads );
	} // for N

	free(x); free(p); free(y);

	exit(0);

	return ;
}

