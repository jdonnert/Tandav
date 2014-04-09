/* OpenMP parallel sorting functions */

void Qsort(void *const pbase, int nElements, size_t size, 
		int (*cmp) (const void *, const void *));
void Qsort_Index(size_t *p, void *const pbase, int nElements, size_t size, 
		int (*cmp) (const void *, const void *));

void test_sort();
