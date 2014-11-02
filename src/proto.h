#ifndef PROTO_H
#define PROTO_H

#include <stdlib.h> 		// system       
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

#include <mpi.h> 			// parallelisation
#include <omp.h>

#include <gsl/gsl_math.h> 	// GNU scientific library

#include "config.h" 		// holds Config #defines

/*
 * Other Code Parameters
 */

#define CHARBUFSIZE 256L 	// Maximum No. of chars in every char buffer
#define NPARTYPE 6L 		// No of particle types
#define MEM_ALIGNMENT 64L	// byte memory alignment
#define PARTALLOCFACTOR 1.1	// Mem overhead for dynamic inbalance
#define MASTER 0			// Global master MPI thread

typedef float Float;		// type of floating point variables in P
typedef uint32_t ID_t;		// type of particle ID
typedef uint64_t intime_t; 	// type of integer time +

#include "macro.h" 			// macro definitions
#include "unit.h" 			// unit functions
#include "constants.h"		// physical constants
#include "aux.h" 			// auxiliary functions 
#include "memory.h"			// memory management
#include "profile.h"		// time measurement
#include "log.h"			// run logging
#include "signal.h"			// signal handlers
#include "sort.h" 			// sort functions
#include "peano.h" 			

/* 
 * Global function prototypes 
 */

extern void Finish();

extern void Print_compile_time_settings();

/* 
 * Add your headers here, #ifdefs go into the .h file 
 */

#include "cosmology.h" 			// cosmology functions 
#include "comov.h"			// Comoving coordinates

/* 
 * Workarounds 
 */

double erand48(unsigned short xsubi[3]);
int posix_memalign(void **memptr, size_t alignment, size_t size);

#endif // PROTO_H
