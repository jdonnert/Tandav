#ifndef SELECT_H
#define SELECT_H

#include "includes.h"
#include "sort.h"

Float Median(const int NData, Float *Data);
Float Select(const int kth, const int NData, Float *Data);
void test_median();

#endif // SELECT_H
