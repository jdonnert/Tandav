#include "io.h"

void IO_Write_Restart_File()
{
	Profile("Restart File");

	char fname[] = { "restartfiles/restart." };

	rprintf("\nWriting restart files %s \n", fname);

	Profile("Restart File");

	return ;
}
