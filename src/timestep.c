#include "globals.h"
#include "timestep.h"

#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)
#define N_INT_BINS (sizeof(intime_t) * CHAR_BIT) // the number of bins in the 
	// integer time line is also the number of bits in an integertime

static int timestep2timebin(const double dt);
static void print_timebins();

static float cosmological_timestep(const int ipart);

struct TimeData Time = { 0 };

/* All active particles get a new step that is smaller than or equal to 
 * the largest active bin. We also update the active particle list 
 * and set the fullstep signal */

void Set_New_Timesteps()
{
	Profile("Timesteps");

	int local_bin_min = N_INT_BINS-1, local_bin_max = 0;

	for (int i = 0; i < NActiveParticles; i++) {

		int ipart = ActiveParticleList[i];
		
		float dt = FLT_MAX;

#ifdef GRAVITY
		float dt_cosmo = cosmological_timestep(ipart); 
#endif
		
		dt = fmin(dt, dt_cosmo);
		
		Assert(dt >= Time.StepMin, "Timestep too small !"
				" Ipart=%d, dt=%g, dt_cosmo=%g", ipart, dt, dt_cosmo);

		P[ipart].TimeBin = Imin(timestep2timebin(dt), Time.MaxActiveBin);
		
		local_bin_min = Imin(local_bin_min, P[ipart].TimeBin);
		local_bin_max = Imax(local_bin_max, P[ipart].TimeBin);
	}
	
	int global_bin_min = N_INT_BINS-1;
	
	MPI_Allreduce(&local_bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);
	
	Time.IntStep = Umin(1ULL << global_bin_min, Time.IntEnd-Time.IntCurrent);

	Time.IntNext += Time.IntStep;
	
	Time.Step = Timebin2Timestep(global_bin_min);

	print_timebins();
	
	Sig.Fullstep = false;

	if (Time.IntCurrent == Time.IntFullStep) {
		
		Sig.Fullstep = true;

		int global_bin_max = 0;

		MPI_Allreduce(&local_bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

		Time.IntFullStep = Umin(Time.IntEnd,  
				Time.IntCurrent + (1ULL << global_bin_max) );

		rprintf("Next full step at t = %g, bin = %d \n\n", 
				Integer2PhysicalTime(Time.IntFullStep), global_bin_max);
	}

	Profile("Timesteps");

	return ;
}

/* The timeline is represented by a 64 bit integer, where an increment of one 
 * corresponds to the whole integration time divided by 2^(N_INT_BINS-1) */

void Setup_Time_Integration()
{
	Time.NextSnap = Time.FirstSnap;

	Time.NSnap = (Time.End - Time.Begin) / Time.BetSnap;

	rprintf("Simulation timeline: \n"
			"   start = %g, end = %g, delta = %g, NSnap = %d \n\n", 
			Time.Begin, Time.End, Time.BetSnap, Time.NSnap);

	Assert(Time.NSnap > 0, "Timeline does not seem to produce any outputs");

	Time.IntBeg = 0;
	Time.IntEnd = 1ULL << (N_INT_BINS - 1);

	Time.IntCurrent = Time.IntBeg;

	Time.StepMax = (Time.End - Time.Begin);
	Time.StepMin =  Time.StepMax / (1ULL << (N_INT_BINS - 1) );

	Time.MaxActiveBin = N_INT_BINS - 1;

	size_t nBytes = Task.NpartTotalMax * sizeof(*ActiveParticleList);
	
	ActiveParticleList = Malloc(nBytes);

	NActiveParticles = Task.NpartTotal;

	for (int i = 0; i < NActiveParticles; i++)
		ActiveParticleList[i] = i;
	
	return ;
}

void Make_Active_Particle_List()
{
	Time.MaxActiveBin = COUNT_TRAILING_ZEROS(Time.IntCurrent); 
	
	int i = 0;

	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) 
		if (P[ipart].TimeBin <= Time.MaxActiveBin)
			ActiveParticleList[i++] = ipart; 

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

double Integer2PhysicalTime(const intime_t IntegerTime)
{
	return Time.Begin + IntegerTime * Time.StepMin;
}

/* Convert a timestep to a power 2 based timebin *
 * via ceil(log2(StepMax/dt)) */

static int timestep2timebin(const double dt)
{
	int exponent;

	frexp(Time.StepMax/dt, &exponent);

	return N_INT_BINS - 1 - exponent - 1;
}

static void print_timebins()
{
	int npart[N_INT_BINS] = { 0 };
	
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++)
		npart[P[ipart].TimeBin]++;

	int npart_global[N_INT_BINS] = { 0 };

	MPI_Reduce(npart, npart_global, N_INT_BINS, MPI_INT, MPI_SUM, Sim.Master, 
			MPI_COMM_WORLD);
	
	int imin = -1, imax = -1;

	for (int i = 0; i < N_INT_BINS; i++)
		if (npart_global[i] != 0 && imin < 0) 
			imin = i;

	for (int i = N_INT_BINS-1; i > -1; i--)
		if (npart_global[i] != 0 && imax < 0)
			imax = i;	

	rprintf("Timesteps at time %g, NActive %d \n"
			"   Bin       nGas           nDM A  dt\n", 
			Time.Current, NActiveParticles);

	for (int i = imax; i > Time.MaxActiveBin; i--)
		rprintf("   %2d    %7d        %7d %s  %g \n", 
			i, 0, npart_global[i], " ", Timebin2Timestep(i));

	for (int i = Time.MaxActiveBin; i >= imin; i--)
		rprintf("   %2d    %7d        %7d %s  %g \n", 
			i, 0, npart_global[i], "X", Timebin2Timestep(i));

	rprintf("\n");

	return ;
}

#undef N_INT_BINS
#undef COUNT_TRAILING_ZEROS

