#include "globals.h"
#include "proto.h"

void Finish()
{
	Finish_Profiler();
	
	Finish_Memory_Management();

	MPI_Finalize();

	return;
}
