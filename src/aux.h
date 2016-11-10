#ifndef AUX_H
#define AUX_H

/* Helper monkeys */

#include "includes.h"

void Reallocate_P_Info(const char *, const char *, int, const int*, size_t*);
void Assert_Info(const char *, const char *, int, int64_t, const char *, ...);
void Warn_Info(const char *, const char *, int, int64_t, const char *, ...);

int32_t imin(const int32_t x, const int32_t y); // breaks naming convention :(
int32_t imax(const int32_t x, const int32_t y);
int64_t Imin(const int64_t x, const int64_t y);
int64_t Imax(const int64_t x, const int64_t y);

uint32_t umin(const uint32_t x, const uint32_t y);
uint32_t umax(const uint32_t x, const uint32_t y);
uint64_t Umin(const uint64_t x, const uint64_t y);
uint64_t Umax(const uint64_t x, const uint64_t y);

Float Sign(const Float x);

void Print_Int_Bits(const __uint128_t val, const int length, const int delta);
void Print_Int_Bits32 (const uint32_t);
void Print_Int_Bits64 (const uint64_t);
void Print_Int_Bits128 (const __uint128_t);

int Fread(void *restrict data, const size_t size, const size_t nWanted,
		FILE *stream);
int Fwrite(void *restrict data, const size_t size, const size_t nWrite,
		FILE *stream);

void Reorder_Array_8(const size_t n, void * restrict p_in, 
														size_t * restrict idx);
void Reorder_Array_4(const size_t n, void * restrict p_in, 
														size_t * restrict idx);
void Reorder_Array_Char(const size_t nBytes, const size_t n, 
								void * restrict p_in, size_t  * restrict idx);

#endif // AUX_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
