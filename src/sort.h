#ifndef SORT_H
#define SORT_H

#include <gsl/gsl_heapsort.h>
#include "includes.h"


/* 
 * OpenMP parallel sorting functions 
 */

void Qsort(const int nThreads, void *const pbase, int nElements, size_t size, 
		int (*cmp) (const void *, const void *));
void Qsort_Index(const int nThreads, size_t *p, void *const pbase, 
		int nElements, size_t size, int (*cmp) (const void *, const void *));

void test_sort();

#endif // SORT_H
