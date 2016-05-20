#ifndef VECTOR_H
#define VECTOR_H

#include "includes.h"
#include "sort.h"
#include "domain.h"

void Find_Vectors();
void Setup_Vectors();

extern struct Vector_Data {
	int * restrict First;
	int * restrict Last;
	peanoKey * restrict Key;
} Vec;

int NVec;

#endif // VECTOR_H
