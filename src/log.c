#include "globals.h"
#include "log.h"

static struct Log_File_Pointers {
	FILE * Profile_Balance;
	FILE * Statistics;
} Log = { NULL };

static struct Statistics {
	double Energy;
	double Momentum;
	double Angular_Momentum[3];
	double CoM[3];
} Stat = { 0 };

static void print_statistics();

/*
 * Write out statistics
 */

void Write_Logs()
{
	Profile_Report_Last(Log.Profile_Balance);
	
	print_statistics(Log.Statistics);

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

	Log.Statistics = fopen(fname, "w");

	Assert(Log.Statistics != NULL, "Can open %s for writing", fname);

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
	fclose(Log.Statistics);
	
	} // omp single nowait

	return ;
}

static void print_statistics(FILE * stream)
{
	#pragma omp single nowait
	{
		double e = 0, p = { 0 }, ang_p[3] = { 0 }, com[3] = { 0 };

		for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
			double mpart = P[ipart].Mass;
			double vpart = ALENGTH3(P[ipart].Vel);
			
			e += mpart * vpart*vpart;
			
			p += mpart * vpart;

			ang_p[0] = mpart * (P[ipart].Pos[1]*P[ipart].Vel[2] 
					- P[ipart].Pos[2]*P[ipart].Vel[1]);
			ang_p[1] = mpart * (P[ipart].Pos[2]*P[ipart].Vel[0] 
					- P[ipart].Pos[0]*P[ipart].Vel[2]);
			ang_p[2] = mpart * (P[ipart].Pos[0]*P[ipart].Vel[1] 
					- P[ipart].Pos[1]*P[ipart].Vel[0]);

			com[0] = mpart * P[ipart].Pos[0];
			com[1] = mpart * P[ipart].Pos[1];
			com[2] = mpart * P[ipart].Pos[2];
		}

		MPI_Reduce(&e, &Stat.Energy, 1, MPI_DOUBLE, MPI_SUM, Sim.Master, 
				MPI_COMM_WORLD);
		
		MPI_Reduce(&p, &Stat.Momentum, 1, MPI_DOUBLE, MPI_SUM, Sim.Master, 
				MPI_COMM_WORLD);

		MPI_Reduce(ang_p, &Stat.Angular_Momentum, 3, MPI_DOUBLE, MPI_SUM, 
				Sim.Master, MPI_COMM_WORLD);

		MPI_Reduce(com, &Stat.CoM, 3, MPI_DOUBLE, MPI_SUM, Sim.Master, 
				MPI_COMM_WORLD);

		if (Task.Is_MPI_Master) 
			fprintf(stream,"%g %g %g %g %g %g %g %g", Stat.Energy, 
					Stat.Momentum, Stat.Angular_Momentum[0], 
					Stat.Angular_Momentum[1], Stat.Angular_Momentum[2], 
					Stat.CoM[0], Stat.CoM[1], Stat.CoM[2] );

	} // omp single nowait

	return ;
}
