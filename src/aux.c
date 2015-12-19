#include "globals.h"
#include "particles.h"

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


/* Error Handling, we use variable arguments to be able
 * to print more informative error messages */

void Assert_Info(const char *func, const char *file, int line,
		int64_t expr, const char *errmsg, ...)
{
	if (expr != 0)
        return;

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

	Assert(nRead == nWanted, "Read %zu bytes, but %zu wanted", nRead, nWanted);

	return nRead;
}

int Fwrite(void *restrict data, const size_t size, const size_t nWrite, 
		FILE *stream)
{
	size_t nWritten = fwrite(data, size, nWrite, stream);

	Assert(nWritten == nWrite, "Wrote %zu bytes, but had %zu",nWrite, nWritten);

	return nWritten;
}


