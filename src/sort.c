/* Simple parallel sorting based on subdiving the list
 * Shamelessly copied from glibc, thereby GPL2 */

#include "globals.h"
/* Byte-wise swap two items of size SIZE. */
#define SWAP(a, b, size)						      \
  do									      \
    {									      \
printf("SWP a %d < b %d \n", *((int8_t*)a), *((int8_t*)b)); fflush(stdout); \
      size_t __size = (size);						      \
      char *__a = (a), *__b = (b);					      \
      do								      \
	{								      \
	  char __tmp = *__a;						      \
	  *__a++ = *__b;						      \
	  *__b++ = __tmp;						      \
	} while (--__size > 0);						      \
    } while (0)

#define MAX_THRESH 1 // change to insertion sort for partitions smaller  

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
  {
    char *lo;
    char *hi;
  } stack_node;

#define STACK_SIZE	(CHAR_BIT * sizeof(stack_node))
#define PUSH(low, high)	((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define	POP(low, high) ((void) (--top, (low = top->lo), (high = top->hi)))
#define	STACK_NOT_EMPTY	(stack < top )

void Qsort (void *const pbase, size_t total_elems, size_t size, 
		int (*cmp) (const void *, const void*))
{
	char *base_ptr = (char *) pbase;

char *pb = (char *) pbase;

  	const size_t max_thresh = MAX_THRESH * size;

printf("TEST \n");

  	if (total_elems == 0) // Avoid lossage with unsigned arithmetic below. 
    	return;

printf("%zu  %d \n", total_elems, MAX_THRESH);

  	if (total_elems > MAX_THRESH) {

      	char *lo = base_ptr; 
      	char *hi = &lo[size * (total_elems - 1)];
 	     	stack_node stack[STACK_SIZE];
      	stack_node *top = stack;
		int n_working = 0;

      	PUSH (lo, hi); // init stack

printf("A \n");
		#pragma omp parallel private(hi, lo) shared(n_working,stack,top)
		while (STACK_NOT_EMPTY || n_working) {

			int worker = 0;

			#pragma omp barrier

			#pragma omp critical
			if (STACK_NOT_EMPTY) { 
				POP(lo,hi);
				
				n_working++;
				worker++;
printf("lo=%p hi=%p \n", lo-pb, hi-pb);
			} 

printf("B \n");
			if (!worker)
				continue;

printf("w %d, nw %d \n", worker, n_working);

        	char *left_ptr;
      		char *right_ptr;
			
	  		/* Select median value from among LO, MID, and HI. Rearrange
	     	 * LO and HI so the three values are sorted. This lowers the
	     	 * probability of picking a pathological pivot value and
	     	 * skips a comparison for both the LEFT_PTR and RIGHT_PTR in
	     	 * the while loops. */

	  		char *mid = lo + size * ((hi - lo) / size >> 1);

printf("0 lo=%p mid=%p  hi=%p %d %d %d\n", 
		lo-pb,mid-pb,  hi-pb, *((int8_t*)lo), *((int8_t*)mid), *((int8_t*)hi));

printf("CC %d \n", (*cmp) ((void *) mid, (void *) lo));
	  		if ( (*cmp) ((void *) mid, (void *) lo) < 0)
	    		SWAP(mid, lo, size);

printf("1 lo=%p mid=%p  hi=%p %d %d %d\n", 
		lo-pb,mid-pb,  hi-pb, *((int8_t*)lo), *((int8_t*)mid), *((int8_t*)hi));
	  		if ( (*cmp) ((void *) hi, (void *) mid) < 0)
	    		SWAP (mid, hi, size);
	  		else
	    		goto jump_over;
	  	
printf("2 lo=%p mid=%p  hi=%p %d %d %d \n", 
		lo-pb,mid-pb,  hi-pb, *((int8_t*)lo), *((int8_t*)mid), *((int8_t*)hi));
			if ( (*cmp) ((void *) mid, (void *) lo) < 0)
	    		SWAP (mid, lo, size);
		
			jump_over:;
printf("3 lo=%p mid=%p  hi=%p %d %d %d\n", 
		lo-pb,mid-pb,  hi-pb, *((int8_t*)lo), *((int8_t*)mid), *((int8_t*)hi));

	  		left_ptr  = lo + size;
	  		right_ptr = hi - size;

printf("C lo=%p mid=%p hi=%p \n", lo-pb , mid-pb, hi-pb);
	  		/* Here's the famous ``collapse the walls'' section of quicksort.
	     	 * Gotta like those tight inner loops!  They are the main reason
	    	 * that this algorithm runs much faster than others. */
	  		do {

printf("left=%d right=%d \n", *(int8_t*)left_ptr, *(int8_t *)right_ptr);
				
				while ((*cmp) ((void *) left_ptr, (void *) mid) < 0)
					left_ptr += size;

	      		while ((*cmp) ((void *) mid, (void *) right_ptr) < 0)
					right_ptr -= size;

	    		if (left_ptr < right_ptr) {
					
					SWAP (left_ptr, right_ptr, size);
printf("mid=%d left=%d right=%d \n", 
		*((int8_t *)(mid)), * ((int8_t*)(left_ptr)), 
		* ((int8_t*)(right_ptr)));
					if (mid == left_ptr)
		    			mid = right_ptr;
		  			else if (mid == right_ptr)
		    			mid = left_ptr;
		  
					left_ptr += size;
		  
					right_ptr -= size;
		
				} else if (left_ptr == right_ptr) {
printf("OUT \n");	
		  			left_ptr += size;
		  			right_ptr -= size;
		  
					break;
				}

for (int j = 0; j < 10; j++) 
printf("%d %d \n", j,((int8_t*)pb)[j]);
			} while (left_ptr <= right_ptr);

printf("left=%p right=%p lo=%p mid=%p hi=%p \n", 
		left_ptr-pb, right_ptr-pb, lo-pb, mid-pb, hi-pb);


			/* Push next iterations to the stack if they are larger than 
			 * the threshold size. Always push the large one first */

        	if ((right_ptr - lo) > (hi - left_ptr)) { // Push larger first 
        		
				if ((size_t) (right_ptr - lo) > max_thresh) {  // push right
					#pragma omp critical
        			PUSH (left_ptr, hi);
printf("PUSH left=%p hi=%p \n", left_ptr-pb, hi-pb);
				}
        		if ((size_t) (hi - left_ptr) > max_thresh) { // push left 
					#pragma omp critical
            		PUSH (lo, right_ptr);
printf("PUSH lo=%p right=%p \n", lo-pb, right_ptr-pb);
				}

			} else {

				if ((size_t) (hi - left_ptr) > max_thresh) { // push left 
					#pragma omp critical
            		PUSH (lo, right_ptr);
printf("PUSH left=%p hi=%p \n", left_ptr-pb, hi-pb);
				}

				if ((size_t) (right_ptr - lo) > max_thresh) {  // push right
					#pragma omp critical
        			PUSH (left_ptr, hi);
printf("PUSH lo=%p right=%p \n", lo-pb, right_ptr-pb);
				}
			}		

			#pragma omp atomic
			n_working--;

printf("NEXT \n");
    	}
	}
printf("FINISH\n");
return ;
  	/* Once the BASE_PTR array is partially sorted by quicksort the rest
     * is completely sorted using insertion sort, since this is efficient
     * for partitions below MAX_THRESH size. BASE_PTR points to the beginning
     * of the array to sort, and END_PTR points at the very last element in
     * the array (*not* one beyond it!). */

  
    char *const end_ptr = &base_ptr[size * (total_elems - 1)];
    char *tmp_ptr = base_ptr;
    char *thresh = min(end_ptr, base_ptr + max_thresh);
    char *run_ptr;

    /* Find smallest element in first threshold and place it at the
       array's beginning.  This is the smallest array element,
       and the operation speeds up insertion sort's inner loop. */

    for (run_ptr = tmp_ptr + size; run_ptr <= thresh; run_ptr += size)
      if ((*cmp) ((void *) run_ptr, (void *) tmp_ptr) < 0)
        tmp_ptr = run_ptr;

    if (tmp_ptr != base_ptr)
      SWAP (tmp_ptr, base_ptr, size);

    /* Insertion sort, running from left-hand-side up to right-hand-side.  */

    run_ptr = base_ptr + size;
    while ((run_ptr += size) <= end_ptr)
      {
	tmp_ptr = run_ptr - size;
	while ((*cmp) ((void *) run_ptr, (void *) tmp_ptr) < 0)
	  tmp_ptr -= size;

	tmp_ptr += size;
        if (tmp_ptr != run_ptr)
          {
            char *trav;

	    trav = run_ptr + size;
	    while (--trav >= run_ptr)
              {
                char c = *trav;
                char *hi, *lo;

                for (hi = lo = trav; (lo -= size) >= tmp_ptr; hi = lo)
                  *hi = *lo;
                *hi = c;
              }
          }
      }
  

  return;
}
