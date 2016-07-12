#ifndef SORT_H
#define SORT_H

#include <gsl/gsl_heapsort.h>
#include "includes.h"

/* 
 * OpenMP sorting functions 
 */

void Qsort(void *data, size_t ndata, size_t size, 
		   int (*cmp) (const void *, const void *));

void Qsort_Index(size_t *perm, void * data, size_t ndata, size_t size, 
		         int (*cmp) (const void *, const void *));

void test_sort();

#endif // SORT_H
Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
