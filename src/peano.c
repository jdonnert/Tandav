#include "globals.h"
#include "proto.h"
#include "peano.h"

/* Here we compute the peano keys and reorder the particles */

//const uint32_t ordering_tab[24][8] = {}, 
//				 orientation_tab[24][8] = {};
const uint32_t ordering_tab[4][4] = 
	{ {0,1,3,2},
	  {0,2,3,1},
	  {3,1,0,2},
	  {3,2,0,1}};
const uint32_t orientation_tab[4][4] = 
	{ {1,0,0,2},
	  {0,1,1,3},
	  {3,2,2,0},
	  {2,3,3,1}};

void print_bits(uint64_t x)
{
	printf("\n");

	for (int i = 0; i < 64; i++) {
		
		uint64_t val = (x & (1 << 63-i)) >> 63-i;
	
		printf("%i", val);

		if (! ((i+1)%4) && i && i!=63)
			printf("_");
	}
		
	printf("\n");

	return ;
}

void print_bits32(uint32_t x)
{
	printf("\n");

	for (int i = 0; i < 32; i++) {
		
		uint64_t val = (x & (1 << 31-i)) >> 31-i;
	
		printf("%i", val);

		if (! ((i+1)%4) && i && i!=31)
			printf("_");
	}
		
	printf("\n");

	return ;
}

void print_bits_peano(uint64_t x)
{
	printf("\n");

	for (uint64_t i = 1; i < 64; i++) {
		
		uint64_t val = (x & (1ULL << 63ULL-i)) >> 63ULL-i;
	
		printf("%i", val);

		if (! ( (i) % 3ULL) && i && i!=63ULL)
			printf("_");
	}
		
	printf("\n");

	return ;
}

/* Based on 'Dynamic Octree Load Balancing Using Space-Filling Curves' 
 * by Campbell et al 2003,
 * Note that the Gadget implementation uses a similar alg.  
 * the old Gadget-2 alg. seems to be from Jagadish 1990 */
peanokey compute_peano_key(const float x, const float y, const float z, 
		const double *boxsize)
{
	uint32_t a = ( x / boxsize[0]) * 0x80000000,
			 b = ( y / boxsize[1]) * 0x40000000,
			 c = ( z / boxsize[2]) * 0x20000000;

	peanokey key = 0;

	//print_bits32(a);
	//print_bits32(b);
	//print_bits32(c);

	const int nBitsPerDim = (int) ( sizeof(peanokey) * CHAR_BIT / 3 ); // 21
	uint32_t or = 0;

	for (int i = 0; i < nBitsPerDim; i++) {

		uint32_t tmp = (a & 0x80000000) | (b & 0x40000000) | (c & 0x20000000);

		tmp >>= 29; // should be one cycle on modern CPUs
		
		//printf("%u %u %u \n", i, tmp, or);

		//key |= tmp;
		//or = orientation_tab[i][0];
		//print_bits_peano(key);

		key |= ordering_tab[or][tmp];
		or = orientation_tab[or][tmp];

		key <<= 3; a <<= 1; b <<= 1; c <<= 1;
	}

	return key;
}



void test()
{
	const double boxsize[] = { 1, 1 , 1};

	const float x = 0, 
		  		y = 0.234234, 
		  		z = 0.892342;

	peanokey key = compute_peano_key( 0, y,z, boxsize);
	print_bits_peano(key);

	printf("\n");

	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++) {
			
			key = compute_peano_key( 0, i*0.25, j*0.25, boxsize);

			printf("%g %g %llu \n", i*0.25, j*0.25, key );

		}

	exit(0);

	return ;
}
