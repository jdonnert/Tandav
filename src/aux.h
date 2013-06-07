inline MPI_Comm Create_MPI_Communicator(const int firstTask, const int lastTask);
void *malloc_info(size_t,const char*,const char*,const int);
void *realloc_info(void *, size_t, const char*, const char*, const int);
void Reallocate_P(size_t*,int);
void safe_free(void *); 
void Assert_Info(const char *, const char *, int, int, const char *, ...);



