/* Simple parallel sorting based on subdiving the list
 * Shamelessly hacked from glibc, thereby GPL2 
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience; 
 * Vol. 23 (11), 1249-1265, 1993. */

#include <gsl/gsl_heapsort.h>
#include "globals.h"

#define PARALLEL_THRESHOLD 10000 // use serial sort below this limit
#define INSERTION_SORT_THRESHOLD 8

#define SWAP(a, b, size)     		\
  	do							    \
    {								\
      size_t __size = (size);		\
      char *__a = (a), *__b = (b);	\
      do							\
	{								\
	  char __tmp = *__a;			\
	  *__a++ = *__b;				\
	  *__b++ = __tmp;				\
	} while (--__size > 0);			\
    } while (0)
	
#define SWAP_SIZE_T(a,b)				\
	do {								\
		size_t *__a = (a), *__b = (b);	\
		size_t __tmp = *__a;			\
		*__a = *__b;					\
		*__b = __tmp;					\
	} while (0)				

#define COMPARE_DATA(a,b,size) ((*cmp) ((void *)(data + *a * size), \
	(void *)(data + *b * size))) // compare with non-permutated data

#define STACK_SIZE (CHAR_BIT * sizeof(size_t))

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
	if (nElements < PARALLEL_THRESHOLD || Sim.NThreads == 1) {

		qsort(pbase, nElements, size, cmp);
	
		return ;
	}

	char *base_ptr = (char *) pbase;
	
	/* initial stack node is just the whole array */
	stack_node stack[STACK_SIZE];

	stack[0].lo = base_ptr; 
	stack[0].hi = &base_ptr[size * (nElements - 1)];

	const int nThreads = 2*(Sim.NThreads/2); // use only 2^N threads

#pragma omp parallel shared(stack) num_threads(nThreads)
	{

	const int tID = Task.ThreadID;

	char *hi, *lo; 

	int i = 0; // number of partitions or iterations 

	for (i = 1; i < nThreads; i *= 2) {

	#pragma omp barrier

		#pragma omp for schedule(static, 1)
		for (int j = 0; j < i; j++) {

			lo = stack[j].lo; // pop from stack
			hi = stack[j].hi;

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

		  	char *left_ptr  = lo + size;
	  		char *right_ptr = hi - size;

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
			stack[2*j].lo = lo;
			stack[2*j].hi = right_ptr;

			stack[2*j + 1].lo = left_ptr;
			stack[2*j + 1].hi = hi;
		} // j
    } // i
	
	#pragma omp barrier

	if (tID <= i) { // serial sort on subpartitions

		size_t chunkSize = (stack[tID].hi-stack[tID].lo)/size + 1;

		qsort(stack[tID].lo, chunkSize, size, cmp);
	} 

	} // omp parallel
	
	return;
}

/* This is an OpenMP parallel external sort.
 * We use the same algorithm as above to divide and conquer.
 * In the first stage we make a partition for every thread.
 * In the second stage a qsort on the subpartition and an insertion sort. 
 * Here we are swaping ONLY the permutation array *perm and
 * are comparing ONLY the data array indexed by *perm */

void Qsort_Index(size_t *perm, void *const data, const int nElements, 
		const size_t datasize, int (*cmp) (const void *, const void *))
{
	if (nElements < PARALLEL_THRESHOLD) {
	
		gsl_heapsort_index(perm, data, nElements, datasize, cmp);
		
		return ;
	}

	for (size_t i = 0; i < nElements; i++ ) 
		perm[i] = i; 

	int nThreads = 2*(Sim.NThreads/2); // N threads has to be power of 2

	idx_stack_node public_stack[nThreads]; // one node per thread  

	public_stack[0].lo = perm;  // initial stack node is just the whole array
	public_stack[0].hi = &perm[nElements - 1];

#pragma omp parallel shared(public_stack)  num_threads(nThreads)
	{

	const int tID = Task.ThreadID;

	size_t *hi, *lo; 
	
	/* First stage: subpartitions */

	for (int i = 1; i < nThreads; i *= 2) { 
 
	#pragma omp barrier

		#pragma omp for schedule(static, 1)
		for (int j = 0; j < i; j++) { // split

			lo = public_stack[j].lo; // pop from public stack
			hi = public_stack[j].hi;

	  		size_t *mid = lo + ((hi - lo) >> 1); // pivot from median

			if (COMPARE_DATA(mid,lo,datasize) < 0) 
				SWAP_SIZE_T (mid, lo);
	
			if (COMPARE_DATA(hi,mid,datasize) < 0)  
				SWAP_SIZE_T (mid, hi);
			else
	    		goto jump_over;
	  	
			if (COMPARE_DATA(mid,lo,datasize) < 0) 
				SWAP_SIZE_T (mid, lo);
			
			jump_over:;
	
			size_t *left  = lo + 1;
	  		size_t *right = hi - 1;

		  	do { // collapse the walls again
			
				while (COMPARE_DATA(left,mid,datasize) < 0) 
					left++;

				while (COMPARE_DATA(mid,right,datasize) < 0)
					right--;

	    		if (left < right) { 
					
					SWAP_SIZE_T (left, right);
					
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

			public_stack[2*j].lo = lo; // Push to public stack
			public_stack[2*j].hi = right;
	
			public_stack[2*j + 1].lo = left;
			public_stack[2*j + 1].hi = hi;
		} // j
    } // i

#pragma omp barrier 

	/* Second stage, every thread: Qsort on subpartitions */

	int top = 0;

	idx_stack_node stack[STACK_SIZE]; // private stack

	stack[top].lo = public_stack[tID].lo;
	stack[top].hi = public_stack[tID].hi;

	const size_t partition_size = stack[top].hi - stack[top].lo + 1;

	if (partition_size < INSERTION_SORT_THRESHOLD)
		goto serial_insertion_sort;

	for (;;) { 

		if (top == -1)
			break;

		lo = stack[top].lo;  // pop from private stack
		hi = stack[top--].hi;

	  	size_t *mid = lo + ((hi - lo) >> 1); 

		if (COMPARE_DATA(mid,lo,datasize) < 0) 
			SWAP_SIZE_T (mid, lo);
	
		if (COMPARE_DATA(hi,mid,datasize) < 0)  
			SWAP_SIZE_T (mid, hi);
		else
	    	goto hop_over;
	  	
		if (COMPARE_DATA(mid,lo,datasize)) 
			SWAP_SIZE_T (mid, lo);
		
		hop_over:;

	  	size_t *left  = lo + 1;
	  	size_t *right = hi - 1;

	  	do { // collapse the walls
		
			while (COMPARE_DATA(left,mid,datasize) < 0) 
				left++;

			while (COMPARE_DATA(mid,right,datasize) < 0)
				right--;

	    	if (left < right) { 
					
				SWAP_SIZE_T (left, right);
					
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

		if (right - lo + 1 > INSERTION_SORT_THRESHOLD) { // push

				stack[++top].lo = lo;
				stack[top].hi = right;
		}
			
		if (hi - left + 1 > INSERTION_SORT_THRESHOLD) {

				stack[++top].lo = left;
				stack[top].hi = hi;
		}
	}
    
	/* Third stage, every thread: insertion sort */

	serial_insertion_sort:;   
	
	size_t *beg = public_stack[tID].lo;
	size_t *end = public_stack[tID].hi;

	size_t *trail = beg;
	size_t *run = NULL;
	
	const size_t *runMax = min(end, beg + INSERTION_SORT_THRESHOLD);

	for (run = beg + 1; run < runMax; run++) // smallest element first
		if (COMPARE_DATA(run, trail, datasize) < 0)
			trail = run;

	if (trail != beg)
		SWAP_SIZE_T(trail, beg);

	run = beg + 1;

	while (++run <= end) { // insertion sort left to right
	
		trail = run - 1;

		while (COMPARE_DATA(run, trail, datasize) < 0)
			trail--;

		trail++;

		if (trail != run)  {

			size_t save = *run;

			hi = run;
			lo = hi - 1;

			while (lo >= trail) 
				*hi-- = *lo--;

			*hi = save;
		}
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
	const size_t N = 500000;
	const size_t Nit = 10;
	int good;

	double *x = malloc(N * sizeof(*x) );
	double *y = malloc(N * sizeof(*y) );
	size_t *p = malloc(N * sizeof(*p) );
	size_t *q = malloc(N * sizeof(*p) );
  	
	double tmp = omp_get_wtime();
  	double tmp2 = omp_get_wtime();

  	double time, time2, time3, deltasum0 = 0, deltasum1 = 0;
	
	for (int j = 0; j < N; j++) 
    	x[j] = y[j] = erand48(Task.Seed);

	for (int i = 0; i < Nit; i++) {
	
		time = omp_get_wtime();

  		Qsort_Index(p, x, N, sizeof(*x), &test_compare);

  		time2 = omp_get_wtime();
	
		gsl_heapsort_index(q, y, N, sizeof(*y), &test_compare);
  		
 	 	time3 = omp_get_wtime();

		deltasum0 += time2-time;
		deltasum1 += time3-time2;
	}

	deltasum0 /= Nit; deltasum1 /= Nit;

  	good = 1;
  
  	for (int i = 1; i < N; i++) {

//		printf("%d %g %g \n", i, x[i], x[p[i]]);

	 	if (x[p[i]] < x[p[i-1]])
			good = 0;
	}

  	if (good)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");
  	
	printf("%zu Parallel:  %g sec, Single:  %g sec, Speedup: %g \n",
		N, deltasum0, deltasum1, deltasum1/deltasum0 );

exit(0);


	/*for (int i = 0; i < Nit; i++) {
	
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

  	good = 1;
  
  	for (int i = 1; i < N; i++) 
	 	if (x[i] < x[i-1])
			good = 0;

  	if (good)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");

  	printf("%zu Parallel:  %g sec, Single:  %g sec, Speedup: %g \n",
		N, deltasum0, deltasum1, deltasum1/deltasum0 );

	exit(0);*/
	
	return ;
}

