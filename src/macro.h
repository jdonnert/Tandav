#if __STDC_VERSION__ < 199901L
# error Recompile with C99 support
#endif

#ifdef DEBUG
#undef MEMORY_MANAGER
#endif

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define MIN(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = MAX(a++, b)
#define MAX(a,b) ((a)>(b)?(a):(b)) // better use functions imin, umin in aux.c

#define ALENGTH3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are slow ! 
#define ALENGTH2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define ASCALPROD3(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are fast !
#define ASCALPROD2(a) (a[0]*a[0] + a[1]*a[1])

#define p2(a) ((a)*(a))
#define p3(a) ((a)*(a)*(a))

#define rprintf(...) if(Task.Is_Master) printf(__VA_ARGS__)
#define mprintf(...) if(Task.Is_MPI_Master) printf(__VA_ARGS__)
#define oprintf(...) if(Task.Is_Thread_Main) printf(__VA_ARGS__)

#define Print_Int_Bits32(x) Print_Int_Bits(x, 32, 3)
#define Print_Int_Bits64(x) Print_Int_Bits(x, 64, 1)
#define Print_Int_Bits128(x) Print_Int_Bits(x, 128, 2)
