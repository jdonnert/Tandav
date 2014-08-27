/* Some useful macros */

#if __STDC_VERSION__ < 199901L
# error Recompile with C99 support
#endif

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0])) 
#define FIELD_SIZEOF(t, f) (sizeof(((t*)0)->f))

#define rprintf(...) if(Task.IsMaster) printf(__VA_ARGS__)

#define min(a,b) ((a)<(b)?(a):(b)) // this doesnt always work: c = max(a++, b)
#define max(a,b) ((a)>(b)?(a):(b))

#define len3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are slow ! 
#define len2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define len3_sq(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define len2_sq(a) (a[0]*a[0] + a[1]*a[1])

#define p2(a) ((a)*(a))  
#define p3(a) ((a)*(a)*(a))
