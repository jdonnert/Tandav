#include "globals.h"

void Finish()
{
	Finish_Profiler();

	Finish_Logs();
	
	Finish_Memory_Management();

	rprintf("\nLive long and prosper\n\n");

	MPI_Finalize();

	return;
}
