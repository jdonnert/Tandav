#include "globals.h"
#include "timestep.h"

#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)

static void make_active_particle_list(const int);
static int timestep2timebin(const double dt);
static void print_timebins();

static float cosmological_timestep(const int ipart);

struct TimeData Time;

static int MaxActiveBin = 63;

/* All active particles get a new step that is smaller than or equal to 
 * the largest active bin. We also update the active particle list 
 * and set the fullstep signal */
void Set_New_Timesteps()
{
	Profile("Timesteps");

	int local_bin_min = 63, local_bin_max = 0;

	#pragma omp parallel for reduction(min:local_bin_min) \
		reduction(max:local_bin_max)
	for (int i = 0; i < NActiveParticles; i++) {
		
		int ipart = ActiveParticleList[i];
		
		float dt = FLT_MAX;

		float dt_cosmo = cosmological_timestep(ipart); 
		
		dt = fmin(dt, dt_cosmo);

		/* Add above */
		Assert(dt > Time.StepMin, "Timestep too small !"
				" Ipart=%d, dt=%g, dt_cosmo=%g", ipart, dt, dt_cosmo);

		int timebin_i = timestep2timebin(dt);

		if (timebin_i < MaxActiveBin)
			P[ipart].TimeBin = timebin_i;

		local_bin_min = Imin(local_bin_min, timebin_i);
		local_bin_max = Imax(local_bin_max, timebin_i);
	}
	
	int global_bin_min = 63;
	
	MPI_Allreduce(&local_bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);

	Time.IntStep = Umin(1ULL << global_bin_min, Time.IntEnd-Time.IntCurrent);
	
	Time.Step = Timebin2Timestep(global_bin_min);

	MaxActiveBin = COUNT_TRAILING_ZEROS(Time.IntCurrent + Time.IntStep); 
	
	make_active_particle_list(MaxActiveBin);

	//print_timebins();

	Flag.Fullstep = false;

	if (Time.IntCurrent == Time.IntSyncPoint) {
		
		Flag.Fullstep = true;

		int global_bin_max = 0;

		MPI_Allreduce(&local_bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);
		
		Time.IntSyncPoint = Umin(Time.IntEnd, 
				Time.IntCurrent + 1ULL << global_bin_max);

		rprintf("Next sync point at %g \n", 
				Integer2PhysicalTime(Time.IntSyncPoint));
	}

	Profile("Timesteps");

	return ;
}

/* The timeline is represented by a 64 bit integer, where an increment of
 * one corresponds to the whole integration time devided by 2^63 */
void Setup_Time_Integration()
{
	Time.NextSnap = Time.FirstSnap;

	Time.NSnap = (Time.End - Time.Begin) / Time.BetSnap;

	rprintf("Simulation timeline: start = %g, end = %g, NSnap = %d \n\n", 
			Time.Begin, Time.End, Time.NSnap);

	Assert(Time.NSnap > 0, "Timeline does not seem to produce any outputs");

	Time.IntBeg = 0x0000000000000000;
	Time.IntEnd = 0x8000000000000000;

	Time.IntCurrent = Time.IntBeg;

	Time.StepMax = Time.End - Time.Begin;
	Time.StepMin =  Time.StepMax / (1ULL << 63);

	size_t nBytes = Sim.NpartTotalMax * sizeof(*ActiveParticleList);
	
	ActiveParticleList = Malloc(nBytes);

	NActiveParticles = Task.NpartTotal;

	for (int i = 0; i < NActiveParticles; i++)
		ActiveParticleList[i] = i;

	return ;
}

static void make_active_particle_list(const int maxActiveBin)
{
	int i = 0;

	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) 
		if (P[ipart].TimeBin <= maxActiveBin)
			ActiveParticleList[i++] = ipart; // all lower bits = 0

	NActiveParticles = i;
	
	return ;
}

/* Cosmological N-body step, Dehnen & Read 2011, eq (21) */
static float cosmological_timestep(const int ipart)
{
	const float acc = len3(P[ipart].Force) / P[ipart].Mass;
	
	float dt = TIME_INT_ACCURACY * sqrt( GRAV_SOFTENING / acc); 
		
	return dt;
}



/* Give the physical timestep from timebin
 * Use Time.Step for the current system step */
double Timebin2Timestep(const int TimeBin)
{
	return Time.StepMin * (1ULL << TimeBin);
}

double Integer2PhysicalTime(const uint64_t IntTime)
{
	return Time.Begin + IntTime * Time.StepMin;
}

/* Convert a timestep to a power 2 based timebin *
 * via ceil(log2(StepMax/dt)) */
static int timestep2timebin(const double dt)
{
	int exponent;

	frexp(Time.StepMax/dt, &exponent);

	return 63 - exponent;
}

static void print_timebins()
{
	int npart[64] = { 0 };
	
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++)
		npart[P[ipart].TimeBin]++;

	int npart_global[64] = { 0 };

	MPI_Reduce(npart, npart_global, 64, MPI_INT, MPI_SUM, MASTER, 
			MPI_COMM_WORLD);
	
	int imin = -1, imax = -1;

	for (int i = 0; i < 64; i++)
		if (npart_global[i] != 0 && imin < 0) 
			imin = i;

	for (int i = 63; i > -1; i--)
		if (npart_global[i] != 0 && imax < 0)
			imax = i;	

	rprintf("Timesteps at time %g \n"
			"   Bin       nGas           nDM   A           dt\n", 
			Time.Current);

	for (int i = imax; i >= imin; i--)
			rprintf("   %2d    %7d        %7d   %d        %g \n", 
					i, 0, npart_global[i], MaxActiveBin/i,Timebin2Timestep(i));

	rprintf("\n");

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}

