#ifdef DEBUG
#undef MEMORY_MANAGER
#endif

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define MIN(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = MAX(a++, b)
#define MAX(a,b) ((a)>(b)?(a):(b)) // better use functions imin, fmin in aux.c

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

#define rprintf(...) if(Task.Is_Master) printf(__VA_ARGS__) // ease printf
#define mprintf(...) if(Task.Is_MPI_Master) printf(__VA_ARGS__)
#define oprintf(...) if(Task.Is_Thread_Main) printf(__VA_ARGS__)

/*
 * Compiler specific stuff
 */

#if __STDC_VERSION__ < 199901L
# error Recompile with C99 support
#endif

#ifdef __INTEL_COMPILER 

#define IVDEP ivdep // IGNORE_VECTOR_DEPENDENCE

#endif // __INTEL_COMPILER

#ifdef _CRAYC

#define IVDEP _CRI ivdep

#endif // _CRAYC
