#include "aux.h"

/*
 * Show bits of an unsigned integer in triplets, 'delta' controls the dot 
 * offset. Has some wrappers associated to print 128, 64, 32 bit peano keys.
 */

void Print_Int_Bits(const __uint128_t val, const int length, const int delta)
{
	__uint128_t one = 1;

	for (int i = length-1; i >= 0; i--) {

		__uint128_t b = val & ( one << i );

		if (b != 0)
			printf("1");
		else
			printf("0");

		if (((i-delta) % 3) == 0 && i != 0)
			printf(".");
	}

	printf("\n");

	fflush(stdout);

	return ;
}

void Print_Int_Bits32 (const uint32_t val)
{
	return Print_Int_Bits((__uint128_t) val, 32, 3);
}
void Print_Int_Bits64 (const uint64_t val)
{
	return Print_Int_Bits((__uint128_t) val, 64, 1);
}
void Print_Int_Bits128 (const __uint128_t val)
{
	return Print_Int_Bits(val, 128, 2);
}

/* 
 * Error Handling, we use variable arguments to print informative messages 
 */

void Assert_Info(const char *func, const char *file, int line,
		int64_t expr, const char *errmsg, ...)
{
	if (expr != 0)
        return;

	Print_Memory_Usage();

	va_list varArgList;

	va_start(varArgList, errmsg);

	fprintf(stderr, "\nERROR (%d:%d) %s:%d : %s() :\n\n	",
			Task.Rank,  Task.Thread_ID, file, line , func);

	vfprintf(stderr, errmsg, varArgList);

	fprintf(stderr, "\n\n");

	fflush(stderr);

	va_end(varArgList);

    MPI_Abort(MPI_COMM_WORLD, -1); // finish him ...

    exit(EXIT_FAILURE); // ... fatality

    return;
}

void Warn_Info(const char *func, const char *file, int line,
		int64_t expr, const char *errmsg, ...)
{
	if (expr == 0)
        return;

	va_list varArgList;

	va_start(varArgList, errmsg);

    fprintf(stderr, "\n\nWARNING (%d:%d): %s:%d : %s() :\n\n  ",
			 Task.Rank, Task.Thread_ID, file, line, func);

	vfprintf(stderr, errmsg, varArgList);

	fprintf(stderr, "\n\n");

	fflush(stderr);

	va_end(varArgList);

    return;
}

/* 
 * Branch free Min/Max functions for signed/unsigned 64/32 bit integers 
 */

int32_t imin(const int32_t x, const int32_t y)
{
  return y ^ ((x ^ y) & -(x < y));
}

int32_t imax(const int32_t x, const int32_t y)
{
  return x ^ ((x ^ y) & -(x < y));
}
uint32_t umin(const uint32_t x, const uint32_t y)
{
  return y ^ ((x ^ y) & -(x < y));
}

uint32_t umax(const uint32_t x, const uint32_t y)
{
  return x ^ ((x ^ y) & -(x < y));
}

int64_t Imin(const int64_t x, const int64_t y)
{
  return y ^ ((x ^ y) & -(x < y));
}

int64_t Imax(const int64_t x, const int64_t y)
{
  return x ^ ((x ^ y) & -(x < y));
}

uint64_t Umin(const uint64_t x, const uint64_t y)
{
  return y ^ ((x ^ y) & -(x < y));
}

uint64_t Umax(const uint64_t x, const uint64_t y)
{
  return x ^ ((x ^ y) & -(x < y));
}

/*
 * Branch free sign function for "Float" type
 */

Float Sign(const Float x)
{
	return ((Float)0 < x) - (x < (Float)0);
}

/*
 * Generic file I/O functions that check for bad reads/writes
 */

int Fread(void *restrict data, const size_t size, const size_t nWanted, 
		FILE *stream)
{
	size_t nRead = fread(data, size, nWanted, stream);

	Assert(nRead == nWanted, "Read %zu objects, but %zu wanted", 
			nRead, nWanted);

	return nRead;
}

int Fwrite(void *restrict data, const size_t size, const size_t nWrite, 
			FILE *stream)
{
	size_t nWritten = fwrite(data, size, nWrite, stream);

	Assert(nWritten == nWrite, "Wrote %zu bytes, but had %zu",
			nWrite, nWritten);

	return nWritten;
}

/*
 * Reorder array p[n] according to idx[n]. idx[] will be changed as well. 
 * We have two versions for 4 and 8 Byte, to avoid using memcpy() and a general
 * implementation using char pointers that takes the element size in bytes.
 */

void Reorder_Array_8(const size_t n, void * restrict p_in, 
									   size_t * restrict idx)
{
	uint64_t * restrict p = (uint64_t * restrict) p_in;

	for (size_t i = 0; i < n; i++) {

   		if (idx[i] == i)
   	    	continue;

		size_t dest = i;

		uint64_t buf = p[i];

		size_t src = idx[i];

 	  	for (;;) {

			p[dest] = p[src];

			idx[dest] = dest;

			dest = src;

			src = idx[dest];

	        if (src == i)
   		        break;
    	}

		p[dest] = buf;

		idx[dest] = dest;
    } // for i

	return ;
}

void Reorder_Array_4(const size_t n, void * restrict p_in, 
					 size_t  * restrict idx)
{
	uint32_t * restrict p = (uint32_t * restrict) p_in;

	for (size_t i = 0; i < n; i++) {

   		if (idx[i] == i)
   	    	continue;

		size_t dest = i;

		uint32_t buf = p[i];

		size_t src = idx[i];

 	  	for (;;) {

			p[dest] = p[src];

			idx[dest] = dest;

			dest = src;

			src = idx[dest];

	        if (src == i)
   		        break;
    	}

		p[dest] = buf;

		idx[dest] = dest;
    } // for i

	return ;
}

void Reorder_Array_Char(const size_t nBytes, const size_t n, 
						void * p_in, size_t  * restrict idx)
{
	char * restrict p = (char * restrict) p_in;

	char buf[nBytes];

	for (size_t i = 0; i < n; i++) {

   		if (idx[i] == i)
   	    	continue;

		size_t dest = i;

		memcpy(buf, p + i*nBytes, nBytes);

		size_t src = idx[i];

 	  	for (;;) {

			memcpy(p + dest*nBytes, p + src*nBytes, nBytes);

			idx[dest] = dest;

			dest = src;

			src = idx[dest];

	        if (src == i)
   		        break;
    	}

		memcpy(p + dest*nBytes, buf, nBytes);

		idx[dest] = dest;
    } // for i

	return ;
}

