#include "globals.h"

void Finish()
{
	Finish_Profiler();

	Finish_Logs();
	
	Finish_Memory_Management();

	rprintf("\nBye\n\n");

	MPI_Finalize();

	return;
}
