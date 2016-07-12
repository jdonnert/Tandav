#include "finish.h"

void Finish()
{
    Finish_Comoving(); // COMOVING

	Finish_Domain_Decomposition();

	Finish_Profiler();

	Finish_Logs();
	
	Finish_Memory_Management();

	rprintf("\nBye\n\n");

	MPI_Finalize();

	return;
}

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
