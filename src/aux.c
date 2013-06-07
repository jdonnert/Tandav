#include "globals.h"


/* Create an MPI Communicator from a range of ranks */
inline MPI_Comm Create_MPI_Communicator(const int firstTask, const int lastTask)
{
	const int mpi_range[1][3] = {readTask, lastTask, 1}; 
	
	MPI_Group mpi_world_group, mpi_read_group;
	
	MPI_Comm mpi_comm_group;
	
	MPI_Comm_group(MPI_COMM_WORLD, &mpi_world_group); 
	
	MPI_Group_range_incl(worldGroup, 1, mpi_range, &mpi_read_group);
	
	MPI_Comm_create(MPI_COMM_WORLD, mpi_read_group, &mpi_comm_read_group);
	
	MPI_Group_free(mpi_world_group);
	
	MPI_Group_free(mpi_read_group);
	
	return mpi_comm_group;
}



/* Memory management */
void *malloc_info(size_t size, const char* file, const char* func, 
		const int line)
{
	void *result = malloc(size);

	Assert_Info(func, file, line, result != NULL || size==0,
			"Allocation failed, %zu Bytes \n" ,size);

	return result;
}

void *realloc_info(void *ptr, size_t size, const char* file, const char* func, 
		const int line)
{
	void * result = realloc(ptr, size);

	Assert_Info(func, file, line, result != NULL || size==0,
			"Reallocation failed: %zu bytes \n" ,size);

	return (result);
}

void safe_free(void *ptr) 
{
    if (ptr != NULL)        
    	free(ptr);
    
    return;
}

/* Reallocates the Particle structures. Takes the relative change
 * as argument, not the total number. Add or Remove via sign argument.
 * Also updates Task.Npart and Task.NPartTotal */
void Reallocate_P(size_t nPart[NO_PART_TYPES], int sign)
{
	for (int type = 0; type < NO_PART_TYPES; type++) {
		Task.Npart[type] += sign * nPart[type];
		Task.NpartTotal += sign * nPart[type];

        Assert(Task.Npart[type]>=0, "Can't allocate negative particles");
	}

    Assert(Task.NpartTotal>=0, "Can't allocate negative particles");

	size_t nBytes = sizeof(*P)*Task.NpartTotal;

	P = realloc(P, nBytes);

	Assert(P!=NULL, "Particle Reallocation failed");

	return;
}

/* Error Handling, we use variable arguments to be able
 * to print more informative error messages */
void Assert_Info(const char *func, const char *file, int line, int expr, 
		const char *errmsg, ...)
{
    if (expr)
        return;

	va_list varArgList;

	va_start(varArgList, errmsg);

	/* we fucked up, tell them */
    fprintf(stderr, "\nTask %d: ERROR %s:%d : %s : \n	", Task.Rank, 
			file, line, func);

	vfprintf(stderr, errmsg, varArgList); 
	
	fprintf(stderr, "\n\n"); 
	
	fflush(stderr);

	va_end(varArgList);

    MPI_Abort(MPI_COMM_WORLD, -1); // finish him ...

    exit(-1); // ... fatality

    return;
}
