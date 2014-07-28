#ifdef MEMORY_MANAGER 
#define Malloc(x) Malloc_info( __func__, __FILE__, __LINE__, x)
#define Realloc(x,y) Realloc_info(__func__, __FILE__, __LINE__, x, y)
#define Free(x) Free_info(__func__, __FILE__, __LINE__, x)
#else
#define Malloc(x) malloc(x)
#define Realloc(x,y) realloc(x,y)
#define Free(x) free(x)
#endif // MEMORY_MANAGER

void *Malloc_info(const char*,const char*,const int, size_t);
void *Realloc_info(const char*, const char*, const int, void *, size_t);
void Free_info(const char* file, const char* func, const int line, void*);

void Init_Memory_Management();
void Print_Memory_Usage();
void Finish_Memory_Management();
void Get_Free_Memory(int *total, int *largest, int *smallest);
