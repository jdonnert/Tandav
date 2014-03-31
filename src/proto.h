#ifndef PROTO_H
#define PROTO_H

#include <stdlib.h> 		// system       
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include <mpi.h> 			// parallelisation
#include <omp.h>

#include <gsl/gsl_math.h> 	// GNU scientific library
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

#include "constants.h"		// physical constants
#include "macro.h" 			// macro definitions
#include "aux.h" 			// auxiliary functions 
#include "unit.h" 			// unit functions
#include "memory.h"			// memory management
#include "profile.h"		// time measurement & logging
#include "sort.h" 			// sort functions

#include "config.h" 		// holds Config #defines

/* Global function prototypes */

extern void Finish();

extern void Print_compile_time_settings();


/* Add your headers here, #ifdefs into .h file */
#include "cosmo.h" 			// cosmology functions 
#include "comov.h"			// Comoving coordinates

/* Workarounds */
double erand48(unsigned short xsubi[3]);
int posix_memalign(void **memptr, size_t alignment, size_t size);
#endif // PROTO_H
