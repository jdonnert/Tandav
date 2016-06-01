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

extern struct Particle_Vector_Blocks{
	int * restrict First;
	int * restrict Last;
} V; // contingouos particle blocks on the same timestep

int NParticle_Vectors;

#endif // VECTOR_H
