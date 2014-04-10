/* Simple parallel sorting based on subdiving the list
 * Shamelessly hacked from glibc, thereby GPL2 
 * Jon Bentley and M. Douglas McIlroy; Software - Practice and Experience; 
 * Vol. 23 (11), 1249-1265, 1993. */

#include "globals.h"
#include <gsl/gsl_heapsort.h>

#define PARALLEL_THRES_QSORT 2000 // use serial sort below this limit
#define PARALLEL_THRES_HEAPS 15000

#define INSERT_THRES 10 // insertion sort threshold
#define SUB_PAR_FAC 8 // sub partitions per thread for load balancing

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

#define STACK_SIZE (CHAR_BIT * sizeof(size_t)) // can't be larger 

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

void Qsort(void *const data_ptr, int nData, size_t size, 
		int (*cmp) (const void *, const void *))
{
	if (nData < PARALLEL_THRES_QSORT || Sim.NThreads == 1) {

		qsort(data_ptr, nData, size, cmp);
	
		return ;
	}

	/* initial stack node is just the whole array */
	stack_node stack[STACK_SIZE] = { NULL };

	stack[0].lo = data_ptr; 
	stack[0].hi = &data_ptr[size * (nData - 1)];

	const int desNumPar = min(Sim.NThreads*SUB_PAR_FAC, nData/INSERT_THRES);

	const size_t minParSize = nData/desNumPar * size;

	int top = 1;

#pragma omp parallel shared(stack,top) 
	{

	for (; top < desNumPar;) {

		int jmax = top;

	#pragma omp barrier

		#pragma omp for 
		for (int j = 0; j < jmax; j++) {

			char *lo = stack[j].lo; // pop from stack
			char *hi = stack[j].hi;

			if (hi - lo < minParSize)  // don't split, too small
				continue;

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
	
			#pragma omp atomic
			top++;

			stack[top - 1].lo = left_ptr;
			stack[top - 1].hi = hi;
		} // j
		
		#pragma omp flush(top)
    } // i
	
	#pragma omp barrier

	#pragma omp for
	for (int i = 0; i < top; i++) { // serial sort on subpartitions

		size_t chunkSize = (stack[i].hi-stack[i].lo)/size + 1;
		
		qsort(stack[i].lo, chunkSize, size, cmp);
	} // i

	} // omp parallel
	
	return;
}

/* This is an OpenMP parallel external sort.
 * We use the same algorithm as above to divide and conquer.
 * In the first stage we sort SUB_PAR_FAC partitions per thread.
 * The second stage a threaded qsort on the subpartitions down,
 * to INSERT_THRES followed by insertion sort on the nearly 
 * ordered subpartition. 
 * Here we are swaping ONLY the permutation array *perm and
 * are comparing ONLY to the data array indexed by *perm */

void Qsort_Index(size_t *perm, void *const data, const int nData, 
		const size_t datasize, int (*cmp) (const void *, const void *))
{
	if (nData < PARALLEL_THRES_HEAPS) { 
	
		gsl_heapsort_index(perm, data, nData, datasize, cmp); // less overhead
		
		return ;
	}

	for (size_t i = 0; i < nData; i++ ) 
		perm[i] = i; 

	const int desNumPar = min(Sim.NThreads, nData/INSERT_THRES);

	const size_t minParSize = nData/desNumPar;
	
	idx_stack_node public_stack[STACK_SIZE] = { NULL };   

	int top = 0;

	public_stack[top].lo = perm;  
	public_stack[top++].hi = &perm[nData - 1];\

#pragma omp parallel shared(public_stack, top)  
	{

	/* First stage: subpartitions, roughly NThreads*SUB_PAR_FAC */

	for (int i = 0; top < desNumPar; i++) { 

		int jmax = top;

	#pragma omp barrier

		#pragma omp for 
		for (int j = 0; j < jmax; j++) { 

			size_t *lo = public_stack[j].lo; // pop from public stack
			size_t *hi = public_stack[j].hi;

			if (hi - lo < minParSize)  // don't split, too small
				continue; 

	  		size_t *mid = lo + ((hi - lo) >> 1); // pivot from median

			if (COMPARE_DATA(mid,lo,datasize) < 0) // sort the three
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

		#pragma omp flush (top)
    } // i

	#pragma omp barrier

	/* Second stage, every thread: Qsort on subpartitions */
	
	#pragma omp for 
	for (int i = 0; i < top; i++) {
	
		idx_stack_node stack[STACK_SIZE] = { NULL }; // private stack

		int next = 0; 

		stack[0].lo = public_stack[i].lo;
		stack[0].hi = public_stack[i].hi;

		const size_t partition_size = stack[0].hi - stack[0].lo + 1;

		if (partition_size < INSERT_THRES)
			goto insertion_sort;

		for (;;) { 

			if (next == -1)
				break;

			size_t *lo = stack[next].lo;  // pop from private stack
			size_t *hi = stack[next--].hi;

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

		insertion_sort:;

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
	
	if (*x > *y)
    	return 1;
	if (*x == *y)
    	return 0;
	else
    	return -1;
}

void test_sort()
{
	const size_t N = 2001;
	const size_t Nit = 200;
	int good;

	double *x = malloc( N * sizeof(*x) );
	double *y = malloc( N * sizeof(*y) );
	size_t *p = malloc( N * sizeof(*p) );
	size_t *q = malloc( N * sizeof(*p) );
  	
	double tmp = omp_get_wtime();
  	double tmp2 = omp_get_wtime();

  	double time, time2, time3, deltasum0 = 0, deltasum1 = 0;
	
	for (int j = 0; j < N; j++) 
    	x[j] = y[j] = erand48(Task.Seed);

/*	for (int i = 0; i < Nit; i++) {
	
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

		//printf("%d %g %g \n", i, x[i], x[p[i]]);

	 	if (x[p[i]] < x[p[i-1]])
			good = 0;
	}

  	if (good)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");
  	
	printf("%zu Parallel:  %g sec, Single:  %g sec, Speedup: %g \n",
		N, deltasum0, deltasum1, deltasum1/deltasum0 );

exit(0);*/


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

	exit(0);

	return ;
}

