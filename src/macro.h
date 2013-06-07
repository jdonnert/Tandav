/* Try to avoid macros (especially global ones). If you can't do without
 * dump them here */

#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

#define rprintf if (!Task.Rank) printf

#define malloc(x) malloc_info(x, __func__, __FILE__, __LINE__)
#define realloc(x,y) realloc_info(x, y, __func__, __FILE__, __LINE__)

#define len3(a) sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) // these are slow ! 
#define len2(a) sqrt(a[0]*a[0] + a[1]*a[1])

#define len3_sq(a) (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
#define len2_sq(a) (a[0]*a[0] + a[1]*a[1])

#define p2(a) (a*a)  
#define p3(a) (a*a*a)

