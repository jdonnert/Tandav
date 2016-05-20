#ifndef INCLUDES_H
#define INCLUDE_H

/*
 * This is the super include file that is present in all other .h files and
 * defines global variables and modules visible everywhere
 */

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

#include <gsl/gsl_math.h>		// GNU scientific library
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "config.h"				// holds Config #defines

#include "globals.h"			// all global variables

/*
 * Global modules - useful everywhere ...
 */

#include "macro.h"				// macro definitions
#include "particles.h"			// particle management    
#include "memory.h"				// memory management
#include "unit.h"				// unit functions
#include "constants.h"			// physical constants
#include "aux.h"				// auxiliary functions 
#include "profile.h"			// time measurement
#include "signal.h"				// signal handlers
#include "vector.h"				// tree leafs as vectors

#endif // INCLUDE_H
