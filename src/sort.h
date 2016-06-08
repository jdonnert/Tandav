#ifndef SORT_H
#define SORT_H

#include <gsl/gsl_heapsort.h>
#include "includes.h"

void Insertion_Sort(void *Data, const size_t nData, const size_t nBytes, 
		int (*cmp) (const void*, const void *));

/* 
 * OpenMP sorting functions 
 */

void Qsort(void *data, size_t nData, size_t size, 
		int (*cmp) (const void *, const void *));

void Qsort_Index(const int nThreads, size_t *p, void *const pbase, 
		int nElements, size_t size, int (*cmp) (const void *, const void *));

void test_sort();

#endif // SORT_H
