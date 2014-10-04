/* Some useful macros */

#if __STDC_VERSION__ < 199901L
# error Recompile with C99 support
#endif

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0])) 
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define rprintf(...) if(Task.Is_Master) printf(__VA_ARGS__) // root print
#define mprintf(...) if(Task.MPI_Rank == MASTER) printf(__VA_ARGS__) // mpi ...

#define MIN(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = max(a++, b)
#define MAX(a,b) ((a)>(b)?(a):(b))

#define ALENGTH3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are slow ! 
#define ALENGTH2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define p2(a) ((a)*(a))  
#define p3(a) ((a)*(a)*(a))
