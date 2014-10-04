#include "globals.h"
#include "log.h"

static struct Log_File_Pointers {
	FILE * Profile_Steps;
} Log = { NULL };

/*
 * Write out statistics
 */

void Write_Logs()
{
	#pragma omp single
	{

	Profile_Report_Last(Log.Profile_Steps);
	
	} // omp master

	return;
}

void Init_Logs()
{
	#pragma omp single
	{
	
	char fname[CHARBUFSIZE] = {""};

	sprintf(fname, "%s/balance", Param.Log_File_Dir);

	Log.Profile_Steps = fopen(fname, "w");

	Assert(Log.Profile_Steps != NULL, "Can open %s for writing", fname);

	} // omp master

	return ;
}


void Finish_Logs()
{
	#pragma omp single
	{
	
	fclose(Log.Profile_Steps);
	
	} // omp master

	return ;
}
