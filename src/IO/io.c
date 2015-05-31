#include "../globals.h"
#include "io.h"

unsigned int Largest_Block_Member_Nbytes()
{
	int imax = 0;

	for (int i = 0; i < NBlocks; i++)
		if (Block[i].Nbytes > Block[imax].Nbytes)
			imax = i;

	return Block[imax].Nbytes;
}

/* 
 * Finds the number of particles for block i 
 */

unsigned int Npart_In_Block(const int i, const int *nPart)
{
	unsigned int list[NPARTYPE + 1] = { 0 };

	for (int i = 0; i < NPARTYPE; i++) {

		list[0] += nPart[i];

		list[i+1] = nPart[i];
	}
	
	int j = (int) Block[i].Target;
	
	return list[j];
}
