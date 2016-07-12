#ifndef MEMORY_H
#define MEMORY_H

#include "includes.h"
#include "timestep.h"
#include "IO/parameter_file.h"

#ifdef MEMORY_MANAGER 
#define Malloc(x,y) Malloc_info(__FILE__, __func__,  __LINE__, x, y)
#define Realloc(x,y,z) Realloc_info(__FILE__, __func__,  __LINE__, x, y, z)
#define Free(x) Free_info( __FILE__, __func__, __LINE__, x)
#define Print_Memory_Usage() Print_Memory_Usage_Info(__FILE__, __func__,  __LINE__)
#else
#define Malloc(x,y) malloc(x)
#define Realloc(x,y,z) realloc(x,y)
#define Free(x) free(x)
#define Print_Memory_Usage() 
#endif // MEMORY_MANAGER

void *Malloc_info(const char*,const char*,const int, size_t, const char*);
void *Realloc_info(const char*, const char*, const int, void *, size_t,
		const char*);
void Free_info(const char* file, const char* func, const int line, void*);

void Init_Memory_Management();
void Print_Memory_Usage_Info(const char* file, const char* func, const int line);
void Finish_Memory_Management();
void Get_Free_Memory(size_t *total, size_t *largest, size_t *smallest);
void *Get_Thread_Safe_Buffer (size_t nBytes);

#endif // MEMORY_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
