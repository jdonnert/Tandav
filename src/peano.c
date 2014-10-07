/* 
 * Here we compute the peano keys and reorder the particles
 */

#include "globals.h"
#include "peano.h"
#include "sort.h"

static peanoKey *keys = NULL;
static size_t *idx = NULL;

int compare_peanokeys(const void * a, const void *b)
{
	const peanoKey *x = (const peanoKey*)a;
	const peanoKey *y = (const peanoKey*)b;

	return (int) (*x - *y);
}

void Sort_Particles_By_Peano_Key()
{
	Profile("Peano-Hilbert order");

	const int npart = Task.Npart_Total;

	#pragma omp single
	{

	if (keys == NULL)
		keys = Malloc(Task.Npart_Total_Max * sizeof(*keys));
	
	if (idx == NULL)
		idx = Malloc(Task.Npart_Total_Max * sizeof(*idx));
	
	} // omp single

	#pragma omp for
	for (int ipart = 0; ipart < npart; ipart++) {

		float x = P[ipart].Pos[0];
		float y = P[ipart].Pos[1];
		float z = P[ipart].Pos[2];

#ifndef PERIODIC
		x += Sim.Boxsize[0] * 0.5;
		y += Sim.Boxsize[1] * 0.5;
		z += Sim.Boxsize[2] * 0.5;
#endif

		keys[ipart] = Peano_Key(x,y,z, Sim.Boxsize);

		P[ipart].Peanokey = keys[ipart];
	}

	Qsort_Index(Sim.NThreads, idx, keys, npart, sizeof(*keys), 
			&compare_peanokeys); 

	#pragma omp for 
	for (int i = 0; i < NActive_Particles; i++)
		Active_Particle_List[i] = idx[Active_Particle_List[i]];

	#pragma omp single
	{

	for (int i = 0; i < npart; i++) {

        if (idx[i] == i)
            continue;

		struct Particle_Data Psrc = P[i];

		int src = i;
		int dest = idx[i];

        for (;;) {

			struct Particle_Data Ptmp = P[dest];
			int idxtmp = idx[dest];

			P[dest] = Psrc;
			idx[dest] = src;

            if (dest == i) 
                break;

			Psrc = Ptmp;
			src = dest = idxtmp;
        }
    }

	} // omp single
	
	Profile("Peano-Hilbert order");

	return ;
}


/* 
 * Construct 64 bit Peano-Hilbert distance in 3D 
 * Yes it's arcane, run as fast as you can.
 * Skilling 2004, AIP 707, 381: "Programming the Hilbert Curve"
 * Note: There is a bug in the code of the paper. See also:
 * Campbell+03 'Dynamic Octree Load Balancing Using Space-Filling Curves' 
 */

peanoKey Peano_Key(const float x, const float y, const float z, 
		const double *boxsize)
{
	const uint32_t m = 0x80000000; // = 1UL << 31;

	uint32_t X[3] = { (y / boxsize[0]) * m, 
				 	  (z / boxsize[1]) * m, 
				      (x / boxsize[2]) * m };

	/* Inverse undo */
    for (uint32_t q = m; q > 1; q >>= 1 ) {

        uint32_t P = q - 1;
        
		if( X[0] & q ) 
			X[0] ^= P;  // invert

        for(int i = 1; i < 3; i++ ) {

			if( X[i] & q ) {

				X[0] ^= P; // invert                              
				
			} else { 
			
				uint32_t t = (X[0] ^ X[i]) & P;  
				
				X[0] ^= t;  
				X[i] ^= t; 
			
			} // exchange
		}
    }

	/* Gray encode (inverse of decode) */
	for(int i = 1; i < 3; i++ )
        X[i] ^= X[i-1];

    uint32_t t = X[2];

    for(int i = 1; i < 32; i <<= 1 )
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for(int i = 1; i >= 0; i-- )
        X[i] ^= t;

	/* bit interleave the transpose array */
	peanoKey key = 0;

	X[1] >>= 1; X[2] >>= 2;	// lowest bits not important

	for (int i = 0; i < 21; i++) {

		uint32_t col = ((X[0] & 0x80000000) 
					  | (X[1] & 0x40000000) 
					  | (X[2] & 0x20000000)) >> 29;
		
		key <<= 3; 

		X[0] <<= 1; X[1] <<= 1; X[2] <<= 1;

		key |= col; 
	} 

	return key;
}

void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--)
		printf("%llu", (val & (1ULL << i) ) >> i);
	
	printf("\n");

	fflush(stdout);

	return ;
}

void test_peanokey()
{
	const double box[3] = { 1.0, 1.0, 1 };
	float a[3] = { 0 };
	int order = 1;
	float delta = 1/pow(2.0, order);
	int n = roundf(1/delta);

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++) 
	for (int k = 0; k < n; k++) {

		a[0] = (i + 0.5) * delta / box[0];
		a[1] = (j + 0.5) * delta / box[1];
		a[2] = (k + 0.5) * delta / box[2];

		peanoKey stdkey =  Peano_Key(a[0], a[1], a[2], box);

		printf("%g %g %g %llu \n", a[0], a[1], a[2], stdkey);

		print_int_bits64(stdkey); printf("\n");

	}
	return ;
}
