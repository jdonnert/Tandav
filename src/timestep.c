#include "globals.h"
#include "timestep.h"

#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)

#define N_INT_BINS (sizeof(intime_t) * CHAR_BIT) // the number of bins in the 
	// integer time line is also the number of bits in intime_t

static int max_active_time_bin();
static void set_particle_timebins(int *bin_max, int *bin_min);
static void set_global_timestep(const int, const int);
static int timestep2timebin(const double dt);
static void make_active_particle_list();
static void print_timebins();

static float cosmological_timestep(const int ipart);

struct TimeData Time = { 0 };

/* 
 * All active particles get a new step that is smaller than or equal to 
 * the largest active bin. We also set the fullstep signal 
 */

void Set_New_Timesteps()
{
	Profile("Timesteps");
#ifndef COMOVING
	rprintf("\nStep <%d> t = %g \n\n", Time.Step_Counter++, Time.Current);
#else
	rprintf("\nStep <%d> a = %g \n\n", Time.Step_Counter++, Time.Current);
#endif
	
	#pragma omp master
	Time.Max_Active_Bin = max_active_time_bin();

	#pragma omp barrier

	int local_bin_min = N_INT_BINS-1; 
	int local_bin_max = 0;

	set_particle_timebins(&local_bin_max, &local_bin_min);

	int global_bin_min = N_INT_BINS-1;
	int global_bin_max = 0;
	
	#pragma omp master
	MPI_Allreduce(&local_bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);

	#pragma omp master
	MPI_Allreduce(&local_bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

	#pragma omp barrier

	set_global_timestep(global_bin_max, global_bin_min);

	make_active_particle_list();

	#pragma omp master
	print_timebins();

	Profile("Timesteps");

	return ;
}

/*
 * The highest active time bin is the last set bit in the current
 * integer time.
 */

static int max_active_time_bin()
{
	if (Int_Time.Current == Int_Time.Beg)
		 return N_INT_BINS-1;

	return COUNT_TRAILING_ZEROS(Int_Time.Current); 
}

/* 
 * Find smallest allowed timestep for local particles 
 * and return local max & min of these bins
 */

static void set_particle_timebins(int *bin_max, int *bin_min)
{
	int local_bin_min = N_INT_BINS-1; 
	int local_bin_max = 0;
	
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		float dt = FLT_MAX;

#ifdef GRAVITY
		float dt_cosmo = cosmological_timestep(ipart); 
		
		dt = fmin(dt, dt_cosmo);
#endif
		// add yours here

		dt = fmin(dt, Time.Step_Max);
		
		Assert(dt >= Time.Step_Min, "Timestep too small or not finite !"
				" Ipart=%d, dt=%g", ipart, dt);

		int want = timestep2timebin(dt);

		int allowed = MAX(Time.Max_Active_Bin, P[ipart].Time_Bin);

		P[ipart].Time_Bin = MIN(want, allowed);

		local_bin_min = MIN(local_bin_min, P[ipart].Time_Bin);

		local_bin_max = MAX(local_bin_max, P[ipart].Time_Bin);
	}

	*bin_min = local_bin_min;
	*bin_max = local_bin_max;

	return ;
}

/*
 * set global timestep and handle fullstep
 */

static void set_global_timestep(const int global_bin_max, 
		const int global_bin_min)
{
	Time.Step = Timebin2Timestep(global_bin_min);
	
	Int_Time.Step = (intime_t) 1 << global_bin_min;
	Int_Time.Step = Umin(Int_Time.Step, Int_Time.End - Int_Time.Current);
	
	Int_Time.Next += Int_Time.Step;
	
	Sig.Fullstep = false;

	if (Int_Time.Current == Int_Time.Full_Step) {
		
		Sig.Fullstep = true;

		Int_Time.Full_Step = Umin(Int_Time.End,  
				Int_Time.Current + ((intime_t) 1 << global_bin_max) );

		rprintf("Next full step at t = %g, bin = %d \n\n", 
				Integer2Physical_Time(Int_Time.Full_Step), global_bin_max);
	}

	return ;
}

/* 
 * The timeline is represented by an integer, where an increment of one 
 * corresponds to the whole integration time divided by 2^(N_INT_BINS-1) 
 */

void Setup_Time_Integration()
{
	Time.Next_Snap = Time.First_Snap;

	Time.NSnap = (Time.End - Time.Begin) / Time.Bet_Snap;

	rprintf("Simulation timeline: \n"
			"   start = %g, end = %g, delta = %g, NSnap = %d \n\n", 
			Time.Begin, Time.End, Time.Bet_Snap, Time.NSnap);

	Assert(Time.NSnap > 0, "Timeline does not seem to produce any outputs");

	Int_Time.Beg = (intime_t) 0;
	Int_Time.End = (intime_t) 1 << (N_INT_BINS - 1);

	Int_Time.Current = Int_Time.Beg;

	Time.Step_Max = (Time.End - Time.Begin);
	Time.Step_Min =  Time.Step_Max / ((intime_t) 1 << (N_INT_BINS - 1) );

	Time.Max_Active_Bin = N_INT_BINS - 1;

	size_t nBytes = Task.Npart_Total_Max * sizeof(*Active_Particle_List);
	
	Active_Particle_List = Malloc(nBytes);

	NActive_Particles = Task.Npart_Total;

	for (int i = 0; i < NActive_Particles; i++) 
		Active_Particle_List[i] = i;

	return ;
}

static void make_active_particle_list()
{
	int i = 0;

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		if (P[ipart].Time_Bin <= Time.Max_Active_Bin)
			Active_Particle_List[i++] = ipart; 
	}

	NActive_Particles = i;

	Assert(NActive_Particles > 0, "No Active Particles, instead %d", i);	

	return ;
}

/* 
 * Give the physical timestep from timebin and vice versa
 * Use Time.Step for the current system step 
 */

double Timebin2Timestep(const int TimeBin)
{
	return Time.Step_Min * ((intime_t) 1 << TimeBin);
}

double Integer2Physical_Time(const intime_t Integer_Time)
{
	return Time.Begin + Integer_Time * Time.Step_Min;
}

/* 
 * Convert a timestep to a power 2 based timebin 
 * via ceil(log2(StepMax/dt)) 
 */

static int timestep2timebin(const double dt)
{
	int exponent;

	frexp(Time.Step_Max/dt, &exponent);

	return N_INT_BINS - 1 - exponent - 1;
}

static void print_timebins()
{
	int npart[N_INT_BINS] = { 0 };
	
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		npart[P[ipart].Time_Bin]++;

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

	rprintf("Systemstep %g, NActive %d \n"
			"   Bin       nGas        nDM A  dt\n", 
			Time.Step, NActive_Particles );

	for (int i = imax; i > Time.Max_Active_Bin; i--)
		rprintf("   %2d    %7d     %7d %s  %0g \n", 
			i, 0, npart_global[i], " ", Timebin2Timestep(i));

	for (int i = Time.Max_Active_Bin; i >= imin; i--)
		rprintf("   %2d    %7d     %7d %s  %0g \n", 
			i, 0, npart_global[i], "X", Timebin2Timestep(i));

	rprintf("\n");

	return ;
}

#undef N_INT_BINS
#undef COUNT_TRAILING_ZEROS

/* 
 * Cosmological N-body step, Dehnen & Read 2011, eq (21) 
 */

static float cosmological_timestep(const int ipart)
{
	const float acc = ALENGTH3(P[ipart].Acc);
	
	float dt = TIME_INT_ACCURACY * sqrt(2*GRAV_SOFTENING / acc); 
		
	return dt;
}
