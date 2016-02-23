#include "globals.h"
#include "timestep.h"
#include "peano.h"
#include "domain.h"
#include "particles.h"

static void reorder_collisionless_particles(const size_t *idx_in);
static void reorder_array_8(const size_t n, void * restrict p_in, 
		size_t  * restrict idx);
static void reorder_array_4(const size_t n, void * restrict p_in, 
		size_t  * restrict idx);

int compare_peanoKeys(const void * a, const void *b)
{
	const peanoKey *x = (const peanoKey *) a;
	const peanoKey *y = (const peanoKey *) b;

	return (int) (*x > *y) - (*x < *y);
}

/* 
 * Here we compute peano Keys and reorder particles
 */

static peanoKey *keys = NULL;
static size_t *idx = NULL;

void Sort_Particles_By_Peano_Key()
{
	Profile("Peano-Hilbert order");
	
	#pragma omp single
	{
	
	keys = Malloc(Task.Npart_Total_Max * sizeof(*keys), "PeanoKeys");

	idx = Malloc(Task.Npart_Total_Max * sizeof(*idx), "Sort Idx");

	} // omp single

	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) 
		keys[ipart] = Peano_Key(P.Pos[0][ipart], P.Pos[1][ipart],
								P.Pos[2][ipart]);

	Qsort_Index(Task.Thread_ID, idx, keys, Task.Npart_Total, sizeof(*keys),
			&compare_peanoKeys);

	#pragma omp single nowait
	Free(keys);

	reorder_collisionless_particles(idx);

	#pragma omp single
 	Free(idx);
	
	Make_Active_Particle_List();
	Make_Active_Particle_Vectors();

	Profile("Peano-Hilbert order");

	return ;
}

static void reorder_collisionless_particles(const size_t *idx_in)
{	
	size_t *idx = NULL;
	size_t nBytes = Task.Npart_Total * sizeof(*idx);

	#pragma omp single
	Assert(Task.Buffer_Size > Task.Npart_Total * sizeof(size_t), 
			"BufferSize >= %g MB for this setup", 
			Sim.NThreads * nBytes/1024.0/1024.0);

	idx = Get_Thread_Safe_Buffer(nBytes);

	#pragma omp for 
	for (int i = 0; i < NP_Fields; i++) { // burn the bus
	
		for (int j = 0; j < P_Fields[i].N; j++) {

			memcpy(idx, idx_in, nBytes); // restore unsorted idx
	
			void * restrict p = Select_Particle(i, j, 0);
	
			if (P_Fields[i].Bytes == 8)
				reorder_array_8(Task.Npart_Total, p, idx);
			else if (P_Fields[i].Bytes == 4)
				reorder_array_4(Task.Npart_Total, p, idx);
			else 
				Assert(false, "Can re-order only 4 or 8 byte entries in P");

		} // for j
	} // for i

	return ;
}


/*
 * Reorder *p according to *idx. *idx will be changed as well. We have two
 * versions for 4 and 8 Byte, to avoid using memcpy().
 */

static void reorder_array_8(const size_t n, void * restrict p_in, 
		size_t  * restrict idx)
{
	uint64_t * restrict p = (uint64_t *) p_in;

	for (size_t i = 0; i < n; i++) {

   		if (idx[i] == i)
   	    	continue;

		size_t dest = i;

		uint64_t buf = p[i];

		size_t src = idx[i];

 	  	for (;;) {

			p[dest] = p[src];

			idx[dest] = dest;

			dest = src;

			src = idx[dest];

	        if (src == i)
   		        break;
    	}

		p[dest] = buf;

		idx[dest] = dest;
    } // for i

	return ;
}

static void reorder_array_4(const size_t n, void * restrict p_in, 
		size_t  * restrict idx)
{
	uint32_t * restrict p = (uint32_t *) p_in;

	for (size_t i = 0; i < n; i++) {

   		if (idx[i] == i)
   	    	continue;

		size_t dest = i;

		uint32_t buf = p[i];

		size_t src = idx[i];

 	  	for (;;) {

			p[dest] = p[src];

			idx[dest] = dest;

			dest = src;

			src = idx[dest];

	        if (src == i)
   		        break;
    	}

		p[dest] = buf;

		idx[dest] = dest;
    } // for i

	return ;
}



/* 
 * Construct a 128 bit Peano-Hilbert distance in 3D, input coordinates 
 * have to be normalized as 0 <= x < 1. Unfortunately it's not clear if 
 * a 128bit type would be portable, though all modern CPUs have 128 bit 
 * registers. Bits 0 & 1 are unused, the most significant bit is 63. Input has
 * to be double for full precision.
 *
 * The normalisation of the position can run into floating point precision 
 * issue for very large domains. Hence we pull the conversion factor 'm' into
 * the conversion first.
 *
 * To be consistent without loss of accuracy all std keys start with the
 * first triplet at level 1, i.e. the last 1 or 2 bits are unused. However the
 * reversed keys carry the zero triplet explicitely, to ease tree traversal.
 *
 * Yes it's totally arcane, run as fast as you can.
 *
 * Skilling 2004, AIP 707, 381: "Programming the Hilbert Curve"
 * Note: There is a bug in the code of the paper. See also:
 * Campbell+03 'Dynamic Octree Load Balancing Using Space-Filling Curves' 
 */

peanoKey Peano_Key(const Float px, const Float py, const Float pz)
{
	const uint64_t m = ((uint64_t) 1) << 63;

	const double fac = m / Domain.Size; // don't run into precision issues ...

	uint64_t X[3] = { (px - Domain.Origin[0]) * fac,
					  (py - Domain.Origin[1]) * fac,
					  (pz - Domain.Origin[2]) * fac };
	/* Inverse undo */

    for (uint64_t q = m; q > 1; q >>= 1 ) {

        uint64_t P = q - 1;

		if( X[0] & q )
			X[0] ^= P;  // invert

        for(int i = 1; i < 3; i++ ) {

			if( X[i] & q ) {

				X[0] ^= P; // invert                              

			} else {

				uint64_t t = (X[0] ^ X[i]) & P;

				X[0] ^= t;
				X[i] ^= t;

			} // exchange
		}
    }

	/* Gray encode (inverse of decode) */

	for(int i = 1; i < 3; i++ )
        X[i] ^= X[i-1];

    uint64_t t = X[2];

    for(int i = 1; i < 64; i <<= 1 )
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for(int i = 1; i >= 0; i-- )
        X[i] ^= t;

	/* branch free bit interleave of transpose array X into key */

	peanoKey key = 0;

	X[1] >>= 1; X[2] >>= 2;	// lowest bits not important

	for (int i = 0; i < N_PEANO_TRIPLETS+1; i++) {

		uint64_t col = ((X[0] & 0x8000000000000000)
					  | (X[1] & 0x4000000000000000)
					  | (X[2] & 0x2000000000000000)) >> 61;

		key <<= 3;
		key |= col;

		X[0] <<= 1;
		X[1] <<= 1;
		X[2] <<= 1;
	}

	key <<= 2;
	
	return key;
}

/*
 * This constructs the peano key with reversed triplet order. The order in the 
 * triplets however is the same ! Also level zero is carried explicitely
 * to ease tree construction. The most significant bits are undefined.
 */

peanoKey Reversed_Peano_Key(const Float px, const Float py, const Float pz)
{
	const uint64_t m = ((uint64_t) 1) << 63;

	const double fac = m / Domain.Size;

	uint64_t X[3] = { (px - Domain.Origin[0]) * fac,
					  (py - Domain.Origin[1]) * fac,
					  (pz - Domain.Origin[2]) * fac };

	/* Inverse undo */

    for (uint64_t q = m; q > 1; q >>= 1) {

        uint64_t P = q - 1;

		if(X[0] & q)
			X[0] ^= P;  // invert

        for(int i = 1; i < 3; i++ ) {

			if(X[i] & q) {

				X[0] ^= P; // invert                              

			} else {

				uint64_t t = (X[0] ^ X[i]) & P;

				X[0] ^= t;
				X[i] ^= t;

			} // exchange
		}
    }

	/* Gray encode (inverse of decode) */

	for(int i = 1; i < 3; i++)
        X[i] ^= X[i-1];

    uint64_t t = X[2];

    for(int i = 1; i < 64; i <<= 1)
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for(int i = 1; i >= 0; i--)
        X[i] ^= t;

	/* branch free reversed (!) bit interleave of transpose array X into key */

	peanoKey key = 0;

	X[0] >>= 18; X[1] >>= 19; X[2] >>= 20;	// lowest bits not important

	for (int i = 0; i < N_PEANO_TRIPLETS+1; i++) {

		uint64_t col = ((X[0] & 0x4) | (X[1] & 0x2) | (X[2] & 0x1));

		key <<= 3;

		key |= col;

		X[0] >>= 1;
		X[1] >>= 1;
		X[2] >>= 1;
	}

	key <<= 3; // include level 0

	return key;
}

/* 
 * Keys of 64 bit length (21 triplets), standard and reversed
 */

shortKey Short_Peano_Key(const Float px, const Float py, const Float pz)
{
	const uint32_t m = ((uint32_t) 1) << 31;

	const double fac = m / Domain.Size;

	uint32_t X[3] = { (px - Domain.Origin[0]) * fac,
					  (py - Domain.Origin[1]) * fac,
					  (pz - Domain.Origin[2]) * fac };

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

	for (int i = 0; i < N_SHORT_TRIPLETS+1; i++) {

		uint64_t col = ((X[0] & 0x80000000)
					  | (X[1] & 0x40000000)
					  | (X[2] & 0x20000000)) >> 29;

		key <<= 3;

		X[0] <<= 1;
		X[1] <<= 1;
		X[2] <<= 1;

		key |= col;
	}

	key <<= 1;

	return key;
}


shortKey Reversed_Short_Peano_Key(const Float px, const Float py, 
		const Float pz)
{
	const uint32_t m = ((uint32_t) 1) << 31;

	const double fac = m / Domain.Size;

	uint32_t X[3] = { (px - Domain.Origin[0]) * fac,
					  (py - Domain.Origin[1]) * fac,
					  (pz - Domain.Origin[2]) * fac };

	/* Inverse undo */

    for (uint32_t q = m; q > 1; q >>= 1) {

        uint32_t P = q - 1;

		if( X[0] & q )
			X[0] ^= P;  // invert

        for(int i = 1; i < 3; i++) {

			if(X[i] & q) {

				X[0] ^= P; // invert                              

			} else {

				uint32_t t = (X[0] ^ X[i]) & P;

				X[0] ^= t;
				X[i] ^= t;

			} // exchange
		}
    }

	/* Gray encode (inverse of decode) */

	for(int i = 1; i < 3; i++)
        X[i] ^= X[i-1];

    uint32_t t = X[2];

    for(int i = 1; i < 32; i <<= 1)
        X[2] ^= X[2] >> i;

    t ^= X[2];

    for(int i = 1; i >= 0; i--)
        X[i] ^= t;

	/* branch free reversed (!) bit interleave of transpose array X into key */

	shortKey key = 0;

	X[1] >>= 1; X[2] >>= 2;	// lowest bits not important

	for (int i = 0; i < N_SHORT_TRIPLETS+1; i++) {

		uint32_t col = ((X[0] & 0x4) | (X[1] & 0x2) | (X[2] & 0x1));

		key <<= 3;

		key |= col;

		X[0] >>= 1;
		X[1] >>= 1;
		X[2] >>= 1;
	}

	key <<= 3; // include level 0

	return key;
}


void Test_Peanokey()
{
	const double box[3]  = { 1, 1, 1};
	Float a[3] = { 0 };
	int order = 2;
	double delta = 1/pow(2.0, order);
	int n = roundf(1/delta);

	a[0] = 0.9999999999; // test one value first
	a[1] = 0.9999999999;
	a[2] = 0.9999999999;

	peanoKey stdkey =  Peano_Key(a[0], a[1], a[2]);
	peanoKey revkey =  Reversed_Peano_Key(a[0], a[1], a[2]);

	printf("a0=%lg a1=%lg a2=%lg 64bit val=%llu  \n", a[0], a[1], a[2], 
			(unsigned long long)(stdkey));

	Print_Int_Bits(stdkey, 128, 2);
	Print_Int_Bits(revkey, 128, 3);

	shortKey stdkey_short = Short_Peano_Key(a[0], a[1], a[2]);
	shortKey revkey_short = Reversed_Short_Peano_Key(a[0], a[1], a[2]);

	Print_Int_Bits(stdkey_short, 64, 1);
	Print_Int_Bits(revkey_short, 64, 3);

	printf("\n");

	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++) 
	for (int k = 0; k < n; k++) {

		a[0] = (i + 0.5) * delta / box[0];
		a[1] = (j + 0.5) * delta / box[1];
		a[2] = (k + 0.5) * delta / box[2];

		peanoKey stdkey =  Peano_Key(a[0], a[1], a[2]);
		peanoKey revkey =  Reversed_Peano_Key(a[0], a[1], a[2]);

		printf("%g %g %g %llu  \n", a[0], a[1], a[2], 
				(unsigned long long)(stdkey));

		Print_Int_Bits128(stdkey);
		Print_Int_Bits(revkey, 128, 3);

		printf("\n");
	}

	exit(0);

	return ;
}

