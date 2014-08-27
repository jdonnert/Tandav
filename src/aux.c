#include "globals.h"

/* branch free Min/Max functions for signed/unsigned 64 bit integers */
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

/* Reallocates the Particle structures. Takes the relative change
 * as argument, not the total number. Add or Remove via sign of nPart.
 * Also updates Task.Npart and Task.NPartTotal. 
 * Expands P so that space for nPart[type] is at offset[type]
 * Contracts P so that the last nPart[type] particles are removed 
 * Note that actually no real allocation is taking place, because
 * that would fragment memory. Instead needs to stay with in the
 * limit set by PARTALLOCFACTOR */
void Reallocate_P_Info(const char *func, const char *file, int line, 
		int dNpart[NPARTYPE], size_t offset_out[NPARTYPE])
{
	if (P == NULL)
		P = Malloc(Sim.NpartTotalMax*sizeof(*P));

	for (int i = 0; i < NPARTYPE; i++)
		Assert(Task.Npart[i] + dNpart[i] <= Sim.NpartMax[i], 
			"Too many particles type %d on this task. \n"
			"Have %d, want %d, max %d \nCurrent PARTALLOCFACTOR = %g", 
			i,Task.Npart[i], dNpart[i], Sim.NpartMax[i], PARTALLOCFACTOR);

	int offset[NPARTYPE] = { 0 }, new_npartTotal = 0;
	int new_npart[NPARTYPE] = { 0 };
	
	for (int type = 0; type < NPARTYPE; type++) { // calc offset

		new_npart[type] = Task.Npart[type] + dNpart[type];
		
		new_npartTotal += new_npart[type];
		
        Assert(new_npart[type] >= 0, "Can't alloc negative particles,"
			" type %d, delta %d, current %d,\n"
			"requested from %s, %s(), line %d", 
			type, dNpart[type], Task.Npart[type], file, func, line);

		if (dNpart[type] == 0) 
			continue; // don't need offset here

		for (int i=0; i <= type; i++) 
			offset[type] += new_npart[i];

		offset[type] -= max(0, dNpart[type]); // correct for dNpart>0
	}
 
	int nMove = Task.NpartTotal; // move particles left
	
	for (int type = 0; type < NPARTYPE; type++) { 

		nMove -= Task.Npart[type];

		if (dNpart[type] >= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type] + fabs(dNpart[type]); 
		int dest = offset[type];

		memmove(&P[dest], &P[src], nMove*sizeof(*P));
	}

	//size_t nBytes = sizeof(*P) * new_npartTotal * PARTALLOCFACTOR;

	nMove = Task.NpartTotal; // move particles right

	for (int type = 0; type < NPARTYPE-1; type++) { 

		nMove -= Task.Npart[type];

		if (dNpart[type] <= 0 || Task.Npart[type] == 0 || nMove == 0)
			continue;

		int src = offset[type];
		int dest = offset[type] + dNpart[type];
		
		memmove(&P[dest], &P[src], nMove*sizeof(*P));
	} 

	Task.NpartTotal = new_npartTotal;
	
	for (int type = 0; type < NPARTYPE; type++) // book keeping
		Task.Npart[type] = new_npart[type];

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

	/* we fucked up, tell them */
    fprintf(stderr, "\nERROR Task %d, Thread %d: file %s, line %d,"
			" %s()\n\n", Task.Rank, Task.ThreadID, file, line, func);

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

	/* Houston, we got a problem here */
    fprintf(stderr, "\nWARNING Task %d, Thread %d: In file %s, function %s(),"
			" line %d :\n\n	", Task.Rank, Task.ThreadID, file, func, line);

	vfprintf(stderr, errmsg, varArgList); 
	
	fprintf(stderr, "\n\n"); 
	
	fflush(stderr);

	va_end(varArgList);

    return;
}


