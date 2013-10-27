#include "../globals.h"
#include "../proto.h"

void Write_Restart_File()
{
	if (Time.Running == Param.TimeLimit) 
		rprintf("Run Time Limit %g min reached.\n", Param.TimeLimit / 60);

	return ;
}
