#include "../globals.h"
#include "../proto.h"
#include "../timestep.h"

void Write_Restart_File()
{
	Profile("Restart File");

	char fname[] = { "restartfiles/restart." };

	rprintf("\nWriting restart files %g \n", fname);

	Profile("Restart File");

	return ;
}
