#ifndef MACRO_H
#define MACRO_H

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0])) // only for arr on stack
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define MIN(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = MAX(a++, b)
#define MAX(a,b) ((a)>(b)?(a):(b)) // better use functions imin, fmin in aux.c

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

#define rprintf(...) if (Task.Is_Master) printf(__VA_ARGS__); fflush(stdout)
#define Profile(x) Profile_Info(__func__, __FILE__, __LINE__, x)
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
#define ALIGN __attribute__((aligned(MEM_ALIGNMENT)))

#endif // __INTEL_COMPILER

#ifdef _CRAYC

#define IVDEP _CRI ivdep
#define ALIGN

#endif // _CRAYC

#endif // MACRO_H
Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
