#include "globals.h"
#include "log.h"
#include "timestep.h"

static struct Log_File_Pointers {
	FILE * Profile_Balance;
	FILE * Properties;
} Log = { NULL };

static void print_statistics();

/*
 * Write out statistics
 */

void Write_Logs()
{
	Profile_Report_Last(Log.Profile_Balance);

	print_statistics(Log.Properties);

	return;
}

void Init_Logs()
{
	if (! Task.Is_MPI_Master)
		return ;

	#pragma omp single nowait
	{

	char fname[CHARBUFSIZE] = {""};

	sprintf(fname, "%s/balance", Param.Log_File_Dir);

	Log.Profile_Balance = fopen(fname, "w");

	Assert(Log.Profile_Balance != NULL, "Can open %s for writing", fname);

	sprintf(fname, "%s/statistics", Param.Log_File_Dir);

	Log.Properties = fopen(fname, "w");

	Assert(Log.Properties != NULL, "Can open %s for writing", fname);

	fprintf(Log.Properties,  
			"# Time Mtot Ekin Momentum[3] Ang_Momentum[3] CoM[3] \n"); 
	
	} // omp single nowait

	return ;
}


void Finish_Logs()
{
	if (! Task.Is_MPI_Master)
		return ;

	#pragma omp single nowait
	{

	fclose(Log.Profile_Balance);
	fclose(Log.Properties);

	} // omp single nowait

	return ;
}

static void print_statistics(FILE * stream)
{
	#pragma omp master
	{
		
	if (Task.Is_MPI_Master)
		fprintf(stream,"%g %g %g %g %g %g %g %g %g %g %g %g \n", Time.Current, 
				Sim.Total_Mass, Sim.Kinetic_Energy, 
				Sim.Momentum[0], Sim.Momentum[1], Sim.Momentum[2], 
				Sim.Angular_Momentum[0], Sim.Angular_Momentum[1], 
				Sim.Angular_Momentum[2],
				Sim.Center_Of_Mass[0], Sim.Center_Of_Mass[1], 
				Sim.Center_Of_Mass[2] );

	} // omp single nowait

	return ;
}
