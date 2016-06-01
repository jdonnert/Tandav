#ifndef MACRO_H
#define MACRO_H

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0])) // works only for non-allocatables
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define MIN(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = MAX(a++, b)
#define MAX(a,b) ((a)>(b)?(a):(b)) // better use functions imin, fmin in aux.c

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

#define Profile(x) Profile_Info(__func__, __FILE__, __LINE__, x)
#define rprintf(...) if (Task.Is_Master) printf(__VA_ARGS__) // root printf
#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Warn(...) Warn_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Reallocate_P(...) \
		Reallocate_P_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

/*
 * Check for some compile time errors
 */

#if __STDC_VERSION__ < 199901L
#error Recompile with C99 support
#endif

#ifdef DEBUG
#undef MEMORY_MANAGER
#endif

/*
 * Compiler specific stuff
 */

#ifdef __INTEL_COMPILER

#define IVDEP ivdep // IGNORE_VECTOR_DEPENDENCE

#endif // __INTEL_COMPILER

#ifdef _CRAYC

#define IVDEP _CRI ivdep

#endif // _CRAYC

#endif // MACRO_H
