#define MAXMEMOBJECTS 99999

void *Malloc_info(const char*,const char*,const int, size_t);
void *Realloc_info(const char*, const char*, const int, void *, size_t);
void Free_info(const char* file, const char* func, const int line,void*);
void Init_Memory_Management();
void Print_Memory_Usage();
void Finish_Memory_Management();
