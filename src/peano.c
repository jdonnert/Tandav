/* 
 * Here we compute the peano Keys and reorder the particles
 */

#include "globals.h"
#include "timestep.h"
#include "peano.h"
#include "domain.h"
#include "sort.h"

#include <gsl/gsl_heapsort.h>

static peanoKey *Keys = NULL;
static size_t *Idx = NULL;

static void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--) {
		printf("%llu", (val & (1ULL << i) ) >> i);
		if (i % 3 == 0 && i != 0)
			printf(".");
	}
	printf("\n");fflush(stdout);

	return ;
}


int compare_peanoKeys(const void * a, const void *b)
{
	const peanoKey *x = (const peanoKey *) a;
	const peanoKey *y = (const peanoKey *) b;

	return (int) (*x > *y) - (*x < *y);
}

static void compute_peano_keys();
static void reorder_collisionless_particles();

void Sort_Particles_By_Peano_Key()
{
	Profile("Peano-Hilbert order");
	
	#pragma omp single
	{

	if (Keys == NULL)
		Keys = Malloc(Task.Npart_Total_Max * sizeof(*Keys), "PeanoKeys");
	
	if (Idx == NULL)
		Idx = Malloc(Task.Npart_Total_Max * sizeof(*Idx), "Sort Idx");
	
	} // omp single

	compute_peano_keys();

	Qsort_Index(Sim.NThreads, Idx, Keys, Task.Npart_Total, sizeof(*Keys), 
			&compare_peanoKeys);

	reorder_collisionless_particles();
	
	Make_Active_Particle_List();

	Profile("Peano-Hilbert order");

	return ;
}

static void compute_peano_keys()
{
	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		float px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		float py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		float pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		P[ipart].Peanokey = Keys[ipart] = Peano_Key(px, py, pz);
	}
	return ;
}

static void reorder_collisionless_particles()
{
	#pragma omp single 
	{ 

	for (size_t i = Task.Npart[0]; i < Task.Npart_Total; i++) {

        if (Idx[i] == i)
            continue;

		size_t dest = i;

		struct Particle_Data Ptmp = P[i];

		size_t src = Idx[i];

        for (;;) {

			P[dest] = P[src];

			Idx[dest] = dest;

			dest = src;

			src = Idx[dest];

            if (src == i) 
                break;
        }

		P[dest] = Ptmp;
		Idx[dest] = dest;

    } // for i
	
	} // omp single

	return ;
}

/* 
 * Construct a 64 bit Peano-Hilbert distance in 3D, input coordinates 
 * have to be normed between 0 < x < 1.
 * Yes it's arcane, run as fast as you can.
 * Skilling 2004, AIP 707, 381: "Programming the Hilbert Curve"
 * Note: There is a bug in the code of the paper. See also:
 * Campbell+03 'Dynamic Octree Load Balancing Using Space-Filling Curves' 
 */

peanoKey Peano_Key(const float x, const float y, const float z)
{
#ifdef DEBUG
	Assert(x >= 0 && x <= 1, "X coordinate of out range for PH key %g", x);
	Assert(y >= 0 && y <= 1, "Y coordinate of out range for PH key %g", y);
	Assert(z >= 0 && z <= 1, "Z coordinate of out range for PH key %g", z);
#endif

	const uint32_t m = 1UL << 31; // = 2^31;

	uint32_t X[3] = { y*m, z*m, x*m };

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

	/* branch free bit interleave of transpose array X into key */

	peanoKey key = 0;

	X[1] >>= 1; X[2] >>= 2;	// lowest bits not important

	for (int i = 0; i < 22; i++) {

		uint32_t col = ((X[0] & 0x80000000) 
					  | (X[1] & 0x40000000) 
					  | (X[2] & 0x20000000)) >> 29;
		
		key <<= 3; 

		X[0] <<= 1; 
		X[1] <<= 1; 
		X[2] <<= 1;

		key |= col; 
	} 

	return key;
}



void test_peanokey()
{
	const double box[3]  = { 1.0, 1, 1};
	float a[3] = { 0 };
	int order = 2;
	float delta = 1/pow(2.0, order);
	int n = roundf(1/delta);

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++) 
	for (int k = 0; k < n; k++) {

		a[0] = (i + 0.5) * delta / box[0];
		a[1] = (j + 0.5) * delta / box[1];
		a[2] = (k + 0.5) * delta / box[2];

		peanoKey stdkey =  Peano_Key(a[0], a[1], a[2]);
		peanoKey gdtkey =  0; //Compute_Peano_Key(a[0], a[1], a[2], box);

		printf("%g %g %g %llu %llu \n", a[0], a[1], a[2], stdkey, gdtkey);

		print_int_bits64(stdkey); 
		print_int_bits64(gdtkey); 
		printf("\n");

	}
	return ;
}
