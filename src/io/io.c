#include "../globals.h"
#include "io.h"

unsigned int largest_block_member_nbytes();
unsigned int npart_in_block(const int, const int *);

unsigned int largest_block_member_nbytes()
{
	int imax = 0;

	for (int i = 0; i < NBlocks; i++)
		if (Block[i].Nbytes > Block[imax].Nbytes)
			imax = i;

	return Block[imax].Nbytes;
} 

/* Finds the number of particles for block i, branch free */
unsigned int npart_in_block(const int i, const int *nPart)
{
	unsigned int result = 0;
	
	for (int j = 0; j < NPARTYPE; j++)
		result += nPart[j] * ((1 << j & Block[i].PartBitMask) >> j);

	return result;
}
