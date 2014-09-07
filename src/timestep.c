#include "globals.h"
#include "timestep.h"

#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)
#define N_INT_BINS (sizeof(intime_t) * CHAR_BIT) // the number of bins in the 
	// integer time line is also the number of bits in intime_t

static int timestep2timebin(const double dt);
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

	int local_bin_min = N_INT_BINS-1, local_bin_max = 0;
	
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		float dt = FLT_MAX;

#ifdef GRAVITY
		float dt_cosmo = cosmological_timestep(ipart); 
#endif
		
		dt = fmin(dt, dt_cosmo);
		
		Assert(dt >= Time.StepMin, "Timestep too small or not finite !"
				" Ipart=%d, dt=%g", ipart, dt);

		int want = timestep2timebin(dt);

		int allowed = Imax(Time.MaxActiveBin, P[ipart].Time_Bin);

		P[ipart].Time_Bin = Imin(want, allowed);

		local_bin_min = Imin(local_bin_min, P[ipart].Time_Bin);
		local_bin_max = Imax(local_bin_max, P[ipart].Time_Bin);
	}

	int global_bin_min = N_INT_BINS-1;
	
	#pragma omp master
	MPI_Allreduce(&local_bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);

	#pragma omp wait
	
	Time.Int_Step = Umin(1ULL << global_bin_min, 
			Time.Int_End - Time.Int_Current);

	Time.Int_Next += Time.Int_Step;
	
	Time.Step = Timebin2Timestep(global_bin_min);

	print_timebins();
	
	Sig.Fullstep = false;

	if (Time.Int_Current == Time.Int_Full_Step) {
		
		Sig.Fullstep = true;

		int global_bin_max = 0;

		#pragma omp master
		MPI_Allreduce(&local_bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

		#pragma omp wait

		Time.Int_Full_Step = Umin(Time.Int_End,  
				Time.Int_Current + (1ULL << global_bin_max) );

		rprintf("Next full step at t = %g, bin = %d \n\n", 
				Integer2Physical_Time(Time.Int_FullStep), global_bin_max);
	}

	Profile("Timesteps");

	return ;
}

/* 
 * The timeline is represented by a 64 bit integer, where an increment of one 
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

	Time.Int_Beg = 0;
	Time.Int_End = 1ULL << (N_INT_BINS - 1);

	Time.Int_Current = Time.Int_Beg;

	Time.Step_Max = (Time.End - Time.Begin);
	Time.Step_Min =  Time.Step_Max / (1ULL << (N_INT_BINS - 1) );

	Time.Max_Active_Bin = N_INT_BINS - 1;

	size_t nBytes = Task.Npart_Total_Max * sizeof(*Active_Particle_List);
	
	Active_Particle_List = Malloc(nBytes);

	NActive_Particles = Task.Npart_Total;

	for (int i = 0; i < NActive_Particles; i++)
		Active_Particle_List[i] = i;
	
	return ;
}

void Make_Active_Particle_List()
{
	Time.Max_Active_Bin = COUNT_TRAILING_ZEROS(Time.Int_Current); 
	
	int i = 0;

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		if (P[ipart].Time_Bin <= Time.Max_Active_Bin)
			Active_Particle_List[i++] = ipart; 
	}

	NActive_Particles = i;

	return ;
}



/* Give the physical timestep from timebin
 * Use Time.Step for the current system step */

double Timebin2Timestep(const int TimeBin)
{
	return Time.Step_Min * (1ULL << TimeBin);
}

double Integer2PhysicalTime(const intime_t Integer_Time)
{
	return Time.Begin + Integer_Time * Time.Step_Min;
}

/* Convert a timestep to a power 2 based timebin *
 * via ceil(log2(StepMax/dt)) */

static int timestep2timebin(const double dt)
{
	int exponent;

	frexp(Time.Step_Max/dt, &exponent);

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
			"   Bin       nGas        nDM A  dt\n", 
			Time.Current, NActiveParticles);

	for (int i = imax; i > Time.MaxActiveBin; i--)
		rprintf("   %2d    %7d     %7d %s  %g \n", 
			i, 0, npart_global[i], " ", Timebin2Timestep(i));

	for (int i = Time.MaxActiveBin; i >= imin; i--)
		rprintf("   %2d    %7d     %7d %s  %g \n", 
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
	const float acc = len3(P[ipart].Acc);
	
	float dt = TIME_INT_ACCURACY * sqrt( 2*GRAV_SOFTENING / acc); 
		
	return dt;
}
