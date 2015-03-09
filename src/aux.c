#include "globals.h"



/*
 * Show bits of an unsigned integer in triplets, 'delta' controls the dot 
 * offset. Has some macros associated to print 128, 64, 32 bit peano keys.
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

		if ((i-delta) % 3 == 0 && i != 0)
			printf(".");
	}

	printf("\n");

	fflush(stdout);

	return ;
}

/* 
 * Reallocate the Particle structures. Takes the relative change
 * as argument, not the total number. Add or Remove via sign of nPart.
 * Also updates Task.Npart and Task.NPartTotal. 
 * Expands P so that space for nPart[type] is at offset[type]
 * Contracts P so that the last nPart[type] particles are removed 
 * Note that actually no real allocation is taking place, because
 * that would fragment memory. Instead this needs to stay with in the
 * limit set by PART_ALLOC_FACTOR. 
 */

void Reallocate_P_Info(const char *func, const char *file, int line,
		int dNpart[NPARTYPE], size_t offset_out[NPARTYPE])
{

	#pragma omp single
	for (int i = 0; i < NPARTYPE; i++)
		Assert(Task.Npart[i] + dNpart[i] <= Task.Npart_Max[i],
			"Too many particles type %d on this task. \n"
			"Have %d, want %d, max %d \nCurrent PARTALLOCFACTOR = %g",
			i, Task.Npart[i], dNpart[i], Task.Npart_Max[i], PART_ALLOC_FACTOR);

	int offset[NPARTYPE] = { 0 }, new_npart_total = 0;
	int new_npart[NPARTYPE] = { 0 };

	#pragma omp single
	for (int type = 0; type < NPARTYPE; type++) { // calc offset

		new_npart[type] = Task.Npart[type] + dNpart[type];

		new_npart_total += new_npart[type];

        Assert(new_npart[type] >= 0, "Can't alloc negative particles,"
			" type %d, delta %d, current %d,\n"
			"requested from %s, %s(), line %d",
			type, dNpart[type], Task.Npart[type], file, func, line);

		if (dNpart[type] == 0)
			continue; // don't need offset here

		for (int i=0; i <= type; i++)
			offset[type] += new_npart[i];

		offset[type] -= MAX(0, dNpart[type]); // correct for dNpart>0
	}

	int nMove = Task.Npart_Total; // move particles left

	#pragma omp single
	for (int type = 0; type < NPARTYPE; type++) {

		nMove -= Task.Npart[type];

		if (dNpart[type] >= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type] + fabs(dNpart[type]);
		int dest = offset[type];

		memmove(&P[dest], &P[src], nMove*sizeof(*P));
	}

	nMove = Task.Npart_Total; // move particles right

	#pragma omp single
	for (int type = 0; type < NPARTYPE-1; type++) {

		nMove -= Task.Npart[type];

		if (dNpart[type] <= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type];
		int dest = offset[type] + dNpart[type];

		memmove(&P[dest], &P[src], nMove*sizeof(*P));
	}

	#pragma omp parallel  // book-keeping
	{

	Task.Npart_Total = new_npart_total;

	for (int type = 0; type < NPARTYPE; type++)
		Task.Npart[type] = new_npart[type];
	}

	if (offset_out != NULL) // return ptrs to freed space
		memcpy(offset_out, offset, NPARTYPE*sizeof(*offset));

	return ;
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

	fprintf(stderr, "\nERROR (%d:%d) %s : %d : %s() :\n\n	",
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

    fprintf(stderr, "\nWARNING (%d:%d): In file %s,\n"
			          "              function %s(), line %d :\n\n		",
			 Task.Rank, Task.Thread_ID, file, func, line);

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
	
	Assert(nWritten == nWrite, "Wrote %zu bytes, but had %zu", 
			nWrite, nWritten);

	return nWritten;
}


