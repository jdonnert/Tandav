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
 * Other Code Parameters
 */

#define CHARBUFSIZE 256L	// Maximum No. of bytes in every char buffer
#define NPARTYPE 6L			// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define MASTER 0L			// Global master MPI task

#ifdef DOUBLE_PRECISION
typedef double Float;		// type of floating point variables in P
#define MPI_MYFLOAT MPI_DOUBLE // corresponding MPI communication type macro
#else
typedef float Float;
#define MPI_MYFLOAT MPI_FLOAT
#endif // DOUBLEPRECISION

#ifdef LONG_IDS
typedef uint64_t ID_t;		// type of particle ID
#else
typedef uint32_t ID_t;		// type of particle ID
#endif

typedef uint32_t intime_t;	// type of integer time 
typedef __uint128_t peanoKey; // long peanokey, 42 triplets / levels
typedef uint64_t shortKey;	// short peanokey, 21 triplets / levels

/* 
 * Global function prototypes 
 */

extern void Finish();
extern void Print_compile_time_settings();

/* 
 * Workarounds 
 */

double erand48(unsigned short xsubi[3]);

/* 
 * global module headers  
 */

#include "macro.h"			// macro definitions
#include "unit.h"			// unit functions
#include "constants.h"		// physical constants
#include "aux.h"			// auxiliary functions 
#include "memory.h"			// memory management
#include "profile.h"		// time measurement
#include "log.h"			// run logging
#include "signal.h"			// signal handlers
#include "sort.h"			// sort functions
#include "select.h"			// find median of an array

/* 
 * Add here, #ifdefs go into the .h file 
 */
#include "cosmology.h"		// cosmology functions 
#include "comov.h"			// Comoving coordinates
#include "periodic.h"		// periodic boundary conditions



#endif // PROTO_H
