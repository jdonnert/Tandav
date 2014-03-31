/* Simple parallel sorting based on subdiving the list
 * Shamelessly copied from glibc, thereby GPL2 */

#include "globals.h"

#define SWAP(a, b, size)		\
  	do {				\
		char tmp[size]; 	\
		memcpy(tmp, a, size); 	\
		memcpy(a, b, size); 	\
		memcpy(b, tmp, size);	\
	} while (0)

#define PARALLEL_THRESHOLD 20000 // use serial sort below this limit
#define STACK_SIZE (Sim.NThreads * sizeof(stack_node))

typedef struct // work queue element that holds partitions
{
    char *lo;
    char *hi;
} stack_node;

void test_sort(); 

/* A stack keeps the partition pointers. We assign one thread per 
 * partition until there is a partition for every thread. 
 * Then we continue serial on every partition. There is little
 * need for synchronisation. However the speedup is suboptimal, because the
 * initial iterations are not parallel */

void Qsort (void *const pbase, int nElements, size_t size, 
		int (*cmp) (const void *, const void*))
{
  	if (nElements <= 1) // don't be silly
    		return;
	
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

	if (nElements > PARALLEL_THRESHOLD) { 
	
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

			/* now put all larger/smaller than the pivot*/
			/* on the right/left */
	  		do { 
				
				while ((*cmp)((void *)left_ptr, (void *)mid)<0)
					left_ptr += size;

				while ((*cmp)((void *)mid,(void *)right_ptr)<0)
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

			stack[2*tID+1].lo = left_ptr;
			stack[2*tID+1].hi = hi;
    		}

		#pragma omp barrier
	
		if (tID < i) { // serial sort on subpartitions
			int chunkSize = (stack[tID].hi-stack[tID].lo)/size+1;

			qsort(stack[tID].lo, chunkSize, size, cmp);
		}
	}  

	} // omp parallel
	
	return;
}

int test_compare(const void * a, const void *b) 
{
	const double *ca = (double*)a;
	const double *cb = (double*)b;

	return (int) (*ca - *cb);
}

void test_sort()
{
	const size_t N = 20000000;

	double *x = malloc(N * sizeof(*x) );
	double *y = malloc(N * sizeof(*y) );

  	for (int i = 0; i < N; i++) {

    		x[i] = rand();

	  	y[i] = x[i];
  	}

  	double time = omp_get_wtime();

  	Qsort(x, N, sizeof(*x), &test_compare);

  	double time2 = omp_get_wtime();
	
 	qsort(y, N, sizeof(*y), &test_compare);
  
  	double time4 = omp_get_wtime();

  	int good = 1;
  
  	for (int i = 1; i < N; i++) 
	 	if (x[i] < x[i-1])
			good = 0;

  	if (good)
	  	printf("Array sorted  :-)\n");
  	else 
	  	printf("Array not sorted :-( \n");

  	printf("%zu Parallel:  %g sec, Single:  %g sec\n",
		N, (time2-time), (time4-time2));

exit(0);
	return ;
}

