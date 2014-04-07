#include "../globals.h"
#include "../proto.h"
#include "../timestep.h"

void Write_Restart_File()
{
	rprintf("Run Time Limit %g min reached.\n", Param.RuntimeLimit / 60);

	return ;
}
