#ifndef PROTO_H
#define PROTO_H

#include <stdlib.h>			// system       
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <mpi.h>			// parallelisation
#include <omp.h>

#include <gsl/gsl_math.h>	// GNU scientific library

#include "config.h"			// holds Config #defines

/*
 * Unexposed Code Parameters
 */

#define CHARBUFSIZE 512L	// Maximum No. of bytes in every char buffer
#define NPARTYPE 6L			// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define MASTER 0L			// Global master MPI task for printing

#ifdef DOUBLE_PRECISION

typedef double Float;		// type of floating point variables in P
#define MPI_MYFLOAT MPI_DOUBLE // corresponding MPI communication type macro
#define SQRT sqrt			// type aware square root function for speed

#else

typedef float Float;
#define MPI_MYFLOAT MPI_FLOAT
#define SQRT sqrtf

#endif // ! DOUBLEPRECISION

#ifdef LONG_IDS

typedef uint64_t ID_t;		// type of particle ID

#else

typedef uint32_t ID_t;		

#endif // ! LONG_IDS

typedef uint32_t intime_t;		// type of integer time 

typedef uint64_t shortKey;		// short peanokey, 64 bit = 21 triplets/levels
typedef __uint128_t peanoKey; 	// long peanokey, 128 bit = 42 triplets/levels

/* 
 * Global prototypes 
 */

enum Start_Parameters {
	READ_IC = 0,
	READ_RESTART = 1,
	READ_SNAP = 2,
	DUMP_PARFILE = 10
};

extern void Finish();
extern void Print_Compile_Time_Settings();

/* 
 * Workaround
 */

double erand48(unsigned short xsubi[3]);

/* 
 * Global module headers  
 */

#include "particles.h"		// particle management    
#include "macro.h"			// macro definitions
#include "unit.h"			// unit functions
#include "constants.h"		// physical constants
#include "aux.h"			// auxiliary functions 
#include "memory.h"			// memory management
#include "profile.h"		// time measurement
#include "log.h"			// run logging
#include "signal.h"			// signal handlers
#include "sort.h"			// sort functions
#include "select.h"			// select n-th element and find median of an array

/* 
 * Add here, #ifdefs go into the .h file 
 */

#include "cosmology.h"		// cosmology functions 
#include "comov.h"			// Comoving coordinates
#include "periodic.h"		// periodic boundary conditions


#endif // PROTO_H
