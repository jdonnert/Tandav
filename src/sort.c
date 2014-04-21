/* Simple parallel sorting based on subdiving the list
 * Shamelessly hacked from glibc, thereby GPL2 
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience; 
 * Vol. 23 (11), 1249-1265, 1993. */

#include "globals.h"
#include <gsl/gsl_heapsort.h>

#define PARALLEL_THRES_QSORT 50000 // use serial sort below this limit
#define PARALLEL_THRES_HEAPS 15000

#define INSERT_THRES 8 // insertion sort threshold

#define SWAP(a, b, size)     		\
  	do {						    \
      size_t __size = (size);		\
      char *__a = (a), *__b = (b);	\
      do {							\
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

#define STACK_SIZE (4*CHAR_BIT * sizeof(size_t))  

typedef struct // work queue element that holds partitions
{
    char *lo;
    char *hi;
} stack_node_char;

typedef struct 
{
    size_t *lo;
    size_t *hi;
} stack_node_size_t;


/* A stack keeps the partition pointers. We assign one thread per 
 * partition until there is a partition for every thread. 
 * Then we continue serial on every partition. There is little
 * need for synchronisation. However the speedup is suboptimal, because the
 * initial iterations are not parallel */

void Qsort(const int nThreads, void *const data_ptr, int nData, size_t size, 
		int (*cmp) (const void *, const void *))
{
	if (nData < PARALLEL_THRES_QSORT || nThreads == 1) {

		qsort(data_ptr, nData, size, cmp);
	
		return ;
	}

	/* partition number and size */
	const int desNumPar = 2 * floor(min(nThreads, nData/INSERT_THRES)/2);
	
	/* initial stack node is just the whole array */
	stack_node_char stack[nThreads*2];

	stack[0].lo = (char *)data_ptr; 
	stack[0].hi = &((char *)data_ptr)[size * (nData - 1)];

#pragma omp parallel shared(stack) num_threads(nThreads)
	{

	int delta = desNumPar;

	for (int i = 0; i < log2(desNumPar); i++) {

		delta >>= 1;

		#pragma omp for schedule (static,1)
		for (int j = 0; j < desNumPar; j+=delta<<1) {

			char *lo = stack[j].lo; // pop from stack
			char *hi = stack[j].hi;

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
			
			stack[j].lo = lo;
			stack[j].hi = right_ptr;
	
			stack[j+delta].lo = left_ptr;
			stack[j+delta].hi = hi;
		} // j
    } // i


	#pragma omp for schedule(static,1)
	for (int i = 0; i < desNumPar; i++) { // serial sort on subpartitions

		size_t chunkSize = (stack[i].hi-stack[i].lo)/size + 1;
		
		qsort(stack[i].lo, chunkSize, size, cmp);
	} // i

	} // omp parallel
	
	return ;
}

/* This is an OpenMP parallel external sort.
 * We use the same algorithm as above to divide and conquer.
 * In the first stage we sort 1-4 partitions per thread.
 * The second stage a threaded qsort on the subpartitions down,
 * to INSERT_THRES followed by insertion sort on the nearly 
 * ordered subpartition. 
 * Here we are swaping ONLY the permutation array *perm and
 * are comparing ONLY to the data array indexed by *perm */

void Qsort_Index(const int nThreads, size_t *perm, void *const data, 
		const int nData, const size_t datasize, 
		int (*cmp) (const void *, const void *))
{
	if (nData < PARALLEL_THRES_HEAPS) { // serial GSL is faster
	
		gsl_heapsort_index(perm, data, nData, datasize, cmp); 
		
		return ;
	}

	for (size_t i = 0; i < nData; i++ ) 
		perm[i] = i; 

	const int desNumPar = min(nThreads, nData/INSERT_THRES);

	const size_t minParSize = INSERT_THRES;
	
	stack_node_size_t public_stack[2 * nThreads];   

	int top = 0;

	public_stack[top].lo = perm;  
	public_stack[top++].hi = &perm[nData - 1];

#pragma omp parallel shared(public_stack, top) num_threads(nThreads)
	{

	/* First stage: subpartitions, roughly NThreads */

	for (int i = 0; top < desNumPar; i++) { 

		int jmax = top;

		#pragma omp for schedule(static,1)
		for (int j = 0; j < jmax; j++) { 

			size_t *lo = public_stack[j].lo; 
			size_t *hi = public_stack[j].hi;

			if (hi - lo < minParSize)  
				continue; 

	  		size_t *mid = lo + ((hi - lo) >> 1); 

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
		 
					break;
				}

			} while (left <= right);

			public_stack[j].lo = lo; // Push to replace old position
			public_stack[j].hi = right;
	
			#pragma omp atomic
			top++;

			public_stack[top - 1].lo = left; //Push to top
			public_stack[top - 1].hi = hi;
		} // j
    } // i

	/* Second stage, every thread: Qsort on subpartitions */
	
	#pragma omp for schedule(static,1)
	for (int i = 0; i < top; i++) {
	
		stack_node_size_t stack[STACK_SIZE] = { { NULL, NULL } }; 

		int next = 0; 

		stack[next].lo = public_stack[i].lo;
		stack[next].hi = public_stack[i].hi;

		size_t partition_size = stack[0].hi - stack[0].lo + 1;

		while (partition_size > INSERT_THRES) { 

			if (next == -1)
				break;

			size_t *lo = stack[next].lo;  // now pop from private stack
			size_t *hi = stack[next--].hi;

	  		size_t *mid = lo + ((hi - lo) >> 1); 
		
			if (COMPARE_DATA(mid,lo,datasize) < 0) 
				SWAP_SIZE_T (mid, lo);
	
			if (COMPARE_DATA(hi,mid,datasize) < 0)  
				SWAP_SIZE_T (mid, hi);
			else
		    	goto hop_over;
	  	
			if (COMPARE_DATA(mid,lo,datasize) < 0) 
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

			if (right - lo + 1 > INSERT_THRES) { // push

					stack[++next].lo = lo;
					stack[next].hi = right;
			}
			
			if (hi - left + 1 > INSERT_THRES) {

					stack[++next].lo = left;
					stack[next].hi = hi;
			}
		} // for()
    
		/* Third stage, every thread: insertion sort */

		size_t *beg = public_stack[i].lo;
		size_t *end = public_stack[i].hi;

		size_t *trail = beg;
		size_t *run = NULL;
	
		const size_t *runMax = min(end, beg + INSERT_THRES);

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

				size_t *hi = run;
				size_t *lo = hi - 1;

				while (lo >= trail) 
					*hi-- = *lo--;

				*hi = save;
			}
		} // while
	} // i

	} // omp parallel
	
	return;
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
	const size_t N = 200000;
	const size_t Nit = 10;
	int good;

	double *x = malloc( N * sizeof(*x) );
	double *y = malloc( N * sizeof(*y) );
	size_t *p = malloc( N * sizeof(*p) );
	size_t *q = malloc( N * sizeof(*q) );

  	clock_t time = clock(), time2 = clock(), time3 = clock();
	double deltasum0 = 0, deltasum1 = 0;
	
	/* external / index sort */

	for (int j = 0; j < N; j++) 
    	x[j] = y[j] = erand48(Task.Seed);

	for (int i = 0; i < Nit; i++) {
	
		time = clock();

  		Qsort_Index(Sim.NThreads, p, x, N, sizeof(*x), &test_compare);

  		time2 = clock();
	
		deltasum0 += time2-time;
	}

	deltasum0 /= Nit; 

	for (int i = 0; i < Nit; i++) {
	
  		time2 = clock();
	
		gsl_heapsort_index(q, y, N, sizeof(*y), &test_compare);
  		
		time3 = clock();

		deltasum1 += time3-time2;
	}

	deltasum1 /= Nit;

  	good = 1;
  
  	for (int i = 1; i < N; i++) {

		//printf("%d %g %g \n", i, x[i], x[p[i]]);

	 	if (x[p[i]] < x[p[i-1]])
			good = 0;
	}

  	if (good == 1)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");
  	
	printf("Index: parallel %g sec; Single %g sec; Speedup: %g \n",
		deltasum0/CLOCKS_PER_SEC, 
		deltasum1/CLOCKS_PER_SEC, 	deltasum1/deltasum0 );

	/* in-place sort */

	for (int j = 0; j < N; j++) 
    	x[j] = y[j] = erand48(Task.Seed);

	for (int i = 0; i < Nit; i++) {

		memcpy(x,y,N*sizeof(*x));
	
		time = clock();

  		Qsort(Sim.NThreads, x, N, sizeof(*x), &test_compare);
  		
		time2 = clock();
	
		deltasum0 += time2-time;
	}

	deltasum0 /= Nit; 

	for (int i = 0; i < Nit; i++) {
		
		memcpy(x,y,N*sizeof(*x));

		time2 = clock();
	
 		qsort(x, N, sizeof(*x), &test_compare);
  	
	 	time3 = clock();

		deltasum1 += time3-time2;
	}

	deltasum1 /= Nit;

  	good = 1;
  
  	for (int i = 1; i < N; i++) 
	 	if (x[i] < x[i-1])
			good = 0;

  	if (good == 1)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");

  	printf("In-place: parallel  %g sec, Single:  %g sec, Speedup: %g \n",
		deltasum0/CLOCKS_PER_SEC, 
		deltasum1/CLOCKS_PER_SEC,deltasum1/deltasum0 );

	exit(0);

	return ;
}

