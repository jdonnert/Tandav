#include "profile.h"
#include  "timestep.h"

#define MAXPROFILEITEMS 64		// Max number of profiling marks

static struct Profiling_Object {
	char Name[CHARBUFSIZE];
	double Tbeg;
	double Tend;
	double ThisLast;	// Last iteration this CPU
	double Total;		// Total time over all iterations / All CPUs
	double Min;			// Min time spend here by a CPU
	double Max;			// Max Time spend here by a CPU
	double Mean;		// Mean Time spend here by all CPUs
	double Imbalance;	// Time wasted waiting for the slowest CPU
} Prof[MAXPROFILEITEMS] = { {"", 0} };

static int NProfObjs = 0;
static double Last_Report_Call = 0;

static inline int find_index_from_name(const char *name);
static double measure_time();


void Init_Profiler()
{
	Last_Report_Call = measure_time();

	return ;
}

void Finish_Profiler()
{
	Profile_Report(stdout);

	return ;
}

/* The profiler works using the unique name given at the first call. Upon the 
 * second call with the same name the measurement unit is then stopped again. */

void Profile_Info(const char* file, const char* func, const int line, 
		const char *name)
{
	#pragma omp single
	{

	const int i = find_index_from_name(name);

	if (i == NProfObjs) { // new item, init

		strncpy(Prof[i].Name, name, CHARBUFSIZE);

		NProfObjs++;
	}

	if (Prof[i].Tbeg != 0) { // stop

		Prof[i].Tend = measure_time();

		Prof[i].ThisLast = Prof[i].Tend - Prof[i].Tbeg;

		Prof[i].Total += Prof[i].ThisLast;

		if (i != 0)
			Prof[i].Tbeg = 0;	

#ifdef DEBUG
		printf("\nDEBUG: (%d:%d) ends %s took %g sec \n", Task.Rank, 
				Task.Thread_ID, name, Prof[i].ThisLast); fflush(stdout);
#endif

	} else { // restart

		Prof[i].Tbeg = measure_time();

		Prof[i].Tend = Prof[i].ThisLast = 0;

#ifdef DEBUG
		printf("\nDEBUG: (%d:%d) starts %s \n", 
				Task.Rank, Task.Thread_ID, name); fflush(stdout);
#endif
	}

	} // omp single

	#pragma omp flush

	return ;
}

void Profile_Report(FILE *stream)
{
	#pragma omp single 
	{

	for (int i = 0; i < NProfObjs; i++) { 

		MPI_Reduce(&Prof[i].Total, &Prof[i].Min, 1, MPI_DOUBLE, 
			MPI_MIN, Master, MPI_COMM_WORLD);

		MPI_Reduce(&Prof[i].Total, &Prof[i].Max, 1, MPI_DOUBLE, 
			MPI_MAX, Master, MPI_COMM_WORLD);

		MPI_Reduce(&Prof[i].Total, &Prof[i].Mean, 1, MPI_DOUBLE, 
			MPI_SUM, Master, MPI_COMM_WORLD);

		Prof[i].Mean /= NRank;

		Prof[i].Imbalance += Prof[i].Max - Prof[i].Min;
	}

	if (! Task.Is_MPI_Master) 
		goto skip; 

	double runtime = Runtime();

	double scale = 1; // sec
	char t_unit[CHARBUFSIZE] = { "sec" };

	if (runtime > 60) { // switch to minutes ?

		scale *= 60;

		runtime /= scale;
	
		sprintf(t_unit,"min ");
	} 

	fprintf(stream, "\nProfiler: All sections, total runtime of %g %s\n"
		"                Name       Total              Max       "
		"Min       Mean       Imbalance\n", runtime, t_unit);

	for (int i = 0; i < NProfObjs; i++ )
		fprintf(stream, 
				"%20s    %8.3f (%4.1f%%)   %8.3f  %8.3f  %8.3f   %8.3f\n",
				Prof[i].Name, Prof[i].Total/scale, 
				Prof[i].Total/scale/runtime*100, Prof[i].Max/scale, 
				Prof[i].Min/scale, Prof[i].Mean/scale, 
				Prof[i].Imbalance/scale);

	

	skip:;

	} // omp single 

	#pragma omp barrier

	return ;
}


void Profile_Report_Last(FILE *stream)
{
	#pragma omp single
	{

	const double now = measure_time();

	double min[MAXPROFILEITEMS] = { 0 },
		   max[MAXPROFILEITEMS] = { 0 },
		   mean[MAXPROFILEITEMS] = { 0 },
		   imbalance[MAXPROFILEITEMS] = { 0 };

	for (int i = 1; i < NProfObjs; i++) {

		MPI_Reduce(&Prof[i].ThisLast, &min[i], 1, MPI_DOUBLE,
			MPI_MIN, Master, MPI_COMM_WORLD);

		MPI_Reduce(&Prof[i].ThisLast, &max[i], 1, MPI_DOUBLE,
			MPI_MAX, Master, MPI_COMM_WORLD);

		MPI_Reduce(&Prof[i].ThisLast, &mean[i], 1, MPI_DOUBLE,
			MPI_SUM, Master, MPI_COMM_WORLD);

		mean[i] /= NRank;

		imbalance[i] = max[i] - min[i];
	}

	if (!Task.Is_MPI_Master)
		goto skip;

	double delta_last = now - Last_Report_Call;

	double scale = 1; // sec

	char fullstep[CHARBUFSIZE] = {" "};

	if (Sig.Sync_Point)
		sprintf(fullstep,", Sync_Point");

	if (delta_last > 60) { // switch to minutes ?

		scale *= 60; // min

	} 
	fprintf(stream, "\nStep %d t=%g\n", Time.Step_Counter, Time.Current);

	//fprintf(stream, "                Name"
	//		"      Imbal        Max       Min      Mean\n");

	for (int i = 0; i < NProfObjs; i++ )
		fprintf(stream, "%20s   %8.3f   %8.3f  %8.3f  %8.3f  \n", Prof[i].Name,
				imbalance[i]/scale, max[i]/scale, min[i]/scale, mean[i]/scale);

	Last_Report_Call = now;

	skip:;

	} // omp single 

	return ;
}

double Runtime()
{
	double now = measure_time();

	return (now - Prof[0].Tbeg) ; // in sec
}

static inline int find_index_from_name(const char *name)
{
	int i = 0;

	for (i = 0; i < NProfObjs; i++)
		if (strncmp(name, Prof[i].Name, CHARBUFSIZE) == 0)
			break;

	return i; // may return i = NProfObjs, i.e. new item
}

static double measure_time()
{
	return MPI_Wtime(); // [s]
}

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
