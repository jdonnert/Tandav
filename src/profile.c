#include "globals.h"
#include "profile.h"

#define MAXPROFILEITEMS 999		// Max number of profiling marks

struct AllTaskStats {
	double This;
	double Min;
	double Max;
	double Mean;
};

static struct Profiling_Object {
	char Name[CHARBUFSIZE];
	double Tbeg;
	double Tend;
	struct SimpleStats{
		double This;
		double All;
	} dTotal; // this & max
	struct AllTaskStats dT;
	struct AllTaskStats Mem;
} Prof[MAXPROFILEITEMS];

static int Last = 0;

int find_index_from_name(const char *name);

void Init_Profiler()
{
	Profile("Whole Program");

	return ;
}

void Finish_Profiler() 
{
	Profile("Whole Program");

	Profile_Report();

	return ;
}

void Profile_Info(const char* file, const char* func, const int line, 
		const char *name)
{
	int i = find_index_from_name(name);

	if (i == Last) { // new item start profiling
		
		strncpy(Prof[i].Name, name, CHARBUFSIZE);

		Last++;

		Prof[i].Tbeg = MPI_Wtime();

	} else { // existing item, stop & reduce
		
		Prof[i].Tend = MPI_Wtime();
	
		Prof[i].dT.This = Prof[i].Tbeg - Prof[i].Tend;

		MPI_Reduce(Prof[i].dT.This, Prof[i].dT.Min, 1, MPI_DOUBLE, MPI_MIN, 0,
				MPI_COMM_WORLD);

		MPI_Reduce(Prof[i].dT.This, Prof[i].dT.Max, 1, MPI_DOUBLE, MPI_MAX, 0,
				MPI_COMM_WORLD);

		MPI_Reduce(Prof[i].dT.This, Prof[i].dT.Mean, 1, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);
		Prof[i].dT.Mean /= Sim.NTask;

		Prof[i].dTotal.This += Prof[i].dT.This;
		Prof[i].dTotal.All += Prof[i].dT.Max;
	}

	return ;
}

void Profile_Report()
{
	return ;
}

void Profile_Report_Last()
{

	return ;
}

void Write_Logs()
{
	return;
}

int find_index_from_name(const char *name)
{
	int i = 0;

	for (i = 0; i<Last; i++)
		if (!strncmp(name, Prof[i].Name, CHARBUFSIZE))
			break;

	return i; // may return i=Last, i.e. new item
}
