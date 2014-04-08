/* Simple parallel sorting based on subdiving the list
 * Shamelessly copied from glibc, thereby GPL2 */

#include <gsl/gsl_heapsort.h>
#include "globals.h"

#define SWAP(a, b, size)		\
  	do {						\
		char tmp[size]; 		\
		memcpy(tmp, a, size); 	\
		memcpy(a, b, size); 	\
		memcpy(b, tmp, size);	\
	} while (0)
	
#define PARALLEL_THRESHOLD 10000 // use serial sort below this limit
#define STACK_SIZE (Sim.NThreads * sizeof(stack_node))
#define IDX_STACK_SIZE (Sim.NThreads * sizeof(stack_node))

typedef struct // work queue element that holds partitions
{
    char *lo;
    char *hi;
} stack_node;

typedef struct 
{
    size_t *lo;
    size_t *hi;
} idx_stack_node;


/* A stack keeps the partition pointers. We assign one thread per 
 * partition until there is a partition for every thread. 
 * Then we continue serial on every partition. There is little
 * need for synchronisation. However the speedup is suboptimal, because the
 * initial iterations are not parallel */

void Qsort(void *const pbase, int nElements, size_t size, 
		int (*cmp) (const void *, const void *))
{
	if (nElements < PARALLEL_THRESHOLD) {

		qsort(pbase, nElements, size, cmp);
	
		return ;
	}

	char *base_ptr = (char *) pbase;
	
	/* initial stack node is just the whole array */
	stack_node stack[STACK_SIZE];

	stack[0].lo = base_ptr; 
	stack[0].hi = &base_ptr[size * (nElements - 1)];

#pragma omp parallel shared(stack) 
	{

	const int tID = Task.ThreadID;
	const int itMax = 2*(Sim.NThreads/2); // use only 2^N threads

	char *hi, *lo; 

	int i = 0; // number of partitions or iterations 

	for (i = 1; i < itMax; i *= 2) {
	
#pragma omp barrier

		if (tID >= i)
			continue; // not enough partitions for these
	
		lo = stack[tID].lo; // pop from stack
		hi = stack[tID].hi;

    	char *left_ptr;
      	char *right_ptr;
			
		/* Find pivot element from median and sort the three. 
		 * That helps to prevent the n^2 worst case */
	  	char *mid = lo + size * ((hi - lo) / size >> 1); 

		if ( (*cmp) ((void *) mid, (void *) lo) < 0) 
	   		SWAP(mid, lo, size);
	
		if ( (*cmp) ((void *) hi, (void *) mid) < 0)
	    	SWAP (mid, hi, size);
	  	else
	    	goto jump_over;
	  	
		if ( (*cmp) ((void *) mid, (void *) lo) < 0)
		 	SWAP (mid, lo, size);
		
		jump_over:;

	  	left_ptr  = lo + size;
	  	right_ptr = hi - size;

		/* now put all larger/smaller than the pivot on the right/left */
	  	do { 
				
			while ((*cmp)((void *)left_ptr, (void *)mid) < 0)
				left_ptr += size;

			while ((*cmp)((void *)mid,(void *)right_ptr) < 0)
				right_ptr -= size;

	    	if (left_ptr < right_ptr) { 
					
				SWAP (left_ptr, right_ptr, size);
					
				if (mid == left_ptr)
		   			mid = right_ptr;
				else if (mid == right_ptr)
		   			mid = left_ptr;
		  
				left_ptr += size;
				right_ptr -= size;
		
			} else if (left_ptr == right_ptr) {
		  			
				left_ptr += size;
				right_ptr -= size;
		 
				break;
			}

		} while (left_ptr <= right_ptr);

		/* Push next iterations / partitions to the stack */
		stack[2*tID].lo = lo;
		stack[2*tID].hi = right_ptr;

		stack[2*tID + 1].lo = left_ptr;
		stack[2*tID + 1].hi = hi;

    }
	
	#pragma omp barrier

/*	if (tID <= i) { // serial sort on subpartitions

		int chunkSize = (stack[tID].hi-stack[tID].lo)/size+1;

		qsort(stack[tID].lo, chunkSize, size, cmp);
	} */

	} // omp parallel
	
	return;
}

/* This is an OpenMP parallel external sort
 * we use the same algorithm as above to divide and conquer, 
 * but then use the GSL heapsort. */
void Sort_Index(size_t *p, void *const datap, const int nElements, 
		const size_t dsize, int (*cmp) (const void *, const void *))
{
	if (nElements < PARALLEL_THRESHOLD) { // fallback

		gsl_heapsort_index(p, datap, nElements, dsize, cmp);

		return ;
	}

	/* Initialise permutation */
	for (size_t i = 0; i < nElements; i++ )
		p[i] = i;

	const size_t size = sizeof(*p);

	char *base_ptr = (char *) datap;
	
	/* initial stack node is just the whole array */
	idx_stack_node stack[IDX_STACK_SIZE];

	stack[0].lo = p; 
	stack[0].hi = &p[nElements - 1];

#pragma omp parallel shared(stack) 
	{

	const int tID = Task.ThreadID;
	const int itMax = 2*(Sim.NThreads/2); // use only 2^N threads

	size_t *hi, *lo; 

	int i = 0; // number of partitions or iterations 

	for (i = 1; i < itMax    +    2; i *= 2) {
	#pragma omp barrier

		if (tID >= i)
			continue; // not enough partitions for these
	
		lo = stack[tID].lo; // pop from stack
		hi = stack[tID].hi;

    	size_t *left_ptr;
      	size_t *right_ptr;
			
		/* Find pivot element from median and sort the three. 
		 * That helps to prevent the n^2 worst case */
	  	size_t *mid = lo + ((hi - lo) >> 1); 

		if ( (*cmp) ((void *) (datap +*mid * dsize), 
					 (void *) (datap + *lo * dsize)) < 0) 
			SWAP (mid, lo, size);
	
		if ( (*cmp) ((void *) (datap + *hi * dsize), 
					 (void *) (datap +*mid * dsize)) < 0) 
			SWAP (mid, hi, size);
		else
	    	goto jump_over;
	  	
		if ( (*cmp) ((void *) (datap +*mid * dsize), 
					 (void *) (datap + *lo * dsize)) < 0) 
			SWAP (mid, lo, size);
		
		jump_over:;

	  	left_ptr  = lo + 1;
	  	right_ptr = hi - 1;

		/* now put all larger/smaller than the pivot on the right/left */
	  	do { 
		
			while ((*cmp)((void *)(datap + *left_ptr * dsize), 
						  (void *)(datap + *mid * dsize)) < 0)
				left_ptr++;

			while ((*cmp)((void *)(datap + *mid * dsize), 
						  (void *)(datap + *right_ptr * dsize)) < 0)
				right_ptr--;

	    	if (left_ptr < right_ptr) { 
					
				SWAP (left_ptr, right_ptr, size);
					
				if (mid == left_ptr)
		   			mid = right_ptr;
				else if (mid == right_ptr)
		   			mid = left_ptr;
		  
				left_ptr++;
				right_ptr--;
		
			} else if (left_ptr == right_ptr) {
		  			
				left_ptr++;
				right_ptr--;
		 
				break;
			}

		} while (left_ptr <= right_ptr);

		/* Push next iterations / partitions to the stack */
		stack[2*tID].lo = lo;
		stack[2*tID].hi = right_ptr;

		stack[2*tID + 1].lo = left_ptr;
		stack[2*tID + 1].hi = hi;

    }
	
	#pragma omp barrier

	if (tID <= i) { // serial sort on subpartitions

		int chunkSize = (stack[tID].hi-stack[tID].lo)/size + 1;
		
	}

	} // omp parallel
	
	return;
}








/* testing */
int test_compare(const void * a, const void *b) 
{
	const double *x = (const double*)a;
	const double *y = (const double*)b;
	
	if (*x > *y)
    	return 1;
	if (*x == *y)
    	return 0;
	else
    	return -1;
}

void test_sort()
{
	const size_t N = 10;
	const size_t Nit = 1;

	double *x = malloc(N * sizeof(*x) );
	double *y = malloc(N * sizeof(*y) );
	size_t *p = malloc(N * sizeof(*p) );
  	
	for (int i = 0; i < N; i++) 
    	x[i] = y[i] = erand48(Task.Seed);
	
	for (int i = 0; i < N; i++) 
		printf("%d %zu %g %g \n",i, p[i], x[i],  x[p[i]]);

	//gsl_sort_index(p, x, 1, N);
	Sort_Index(p, x, N, sizeof(*x), &test_compare);
	Qsort(y, N, sizeof(*y), &test_compare);	
	for (int i = 0; i < N; i++) 
		printf("%d %zu %g %g %g \n",i, p[i], x[i],  x[p[i]], y[i]);
	
	exit(0);
	
	double tmp = omp_get_wtime();
  	double tmp2 = omp_get_wtime();

  	double time, time2, time3, deltasum0 = 0, deltasum1 = 0;
	
	for (int i = 0; i < Nit; i++) {
	
		for (int j = 0; j < N; j++) 
    		x[j] = y[j] = erand48(Task.Seed);

		time = omp_get_wtime();

  		Qsort(x, N, sizeof(*x), &test_compare);

  		time2 = omp_get_wtime();
	
 		qsort(y, N, sizeof(*y), &test_compare);
  	
 	 	time3 = omp_get_wtime();

		deltasum0 += time2-time;
		deltasum1 += time3-time2;
	}

	deltasum0 /= Nit; deltasum1 /= Nit;

  	int good = 1;
  
  	for (int i = 1; i < N; i++) 
	 	if (x[i] < x[i-1])
			good = 0;

  	if (good)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");

  	printf("%zu Parallel:  %g sec, Single:  %g sec, Speedup: %g \n",
		N, deltasum0, deltasum1, deltasum1/deltasum0 );

	exit(0);
	
	return ;
}

