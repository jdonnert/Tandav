#ifdef MEMORY_MANAGER 
#define Malloc(x,y) Malloc_info(__FILE__, __func__,  __LINE__, x, y)
#define Realloc(x,y,z) Realloc_info(__FILE__, __func__,  __LINE__, x, y, z)
#define Free(x) Free_info( __FILE__, __func__, __LINE__, x)
#else
#define Malloc(x,y) malloc(x)
#define Realloc(x,y,z) realloc(x,y)
#define Free(x) free(x)
#endif // MEMORY_MANAGER

void *Malloc_info(const char*,const char*,const int, size_t, const char*);
void *Realloc_info(const char*, const char*, const int, void *, size_t,
		const char*);
void Free_info(const char* file, const char* func, const int line, void*);

void Init_Memory_Management();
void Print_Memory_Usage();
void Finish_Memory_Management();
void Get_Free_Memory(int *total, int *largest, int *smallest);
void *Get_Thread_Safe_Buffer (size_t nBytes);
