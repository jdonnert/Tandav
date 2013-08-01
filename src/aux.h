void *Malloc_info(const char*,const char*,const int, size_t);
void *Realloc_info(const char*, const char*, const int, void *, size_t);
void Free_info(const char* file, const char* func, const int line,void*);
void Reallocate_P_Info(const char *, const char *, int, int*, int*);
void safe_free(void *); 
void Assert_Info(const char *, const char *, int, int, const char *, ...);

