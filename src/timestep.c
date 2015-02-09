#include "globals.h"
#include "timestep.h"

static int max_active_time_bin();
static void set_particle_timebins();
static void set_global_timestep();
static int timestep2timebin(const double dt);
static void print_timebins();
static float cosmological_timestep(const int ipart);

/* 
 * The number of bins is given by the number of bits in an integer time 
 */

#define N_INT_BINS (sizeof(intime_t) * CHAR_BIT) 
#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)

struct TimeData Time = { 0 };
struct IntegerTimeLine Int_Time = { 0 };

static int time_bin_min = N_INT_BINS-1, time_bin_max = 0;

/* 
 * All active particles get a new step that is smaller than or equal to 
 * the largest active bin. We also set the fullstep signal 
 */

void Set_New_Timesteps()
{
	Profile("Timesteps");
	
	set_particle_timebins();
	
	#pragma omp single 
	{
	
	MPI_Allreduce(MPI_IN_PLACE, &time_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, &time_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

	} // omp single

	set_global_timestep(time_bin_max, time_bin_min);

	#pragma omp single 
	Time.Max_Active_Bin = max_active_time_bin();

	#pragma omp flush (Time)

	Make_Active_Particle_List();

#ifndef COMOVING
	rprintf("\nStep <%d> t = %g -> %g\n\n", 
			Time.Step_Counter++, Time.Current, 
			Integer2Physical_Time(Int_Time.Next) );
#else
	rprintf("\nStep <%d> a = %g -> %g\n\n", 
			Time.Step_Counter++, Time.Current, Int_Time.Current, 
			Integer2Physical_Time(Int_Time.Next));
#endif

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
	return COUNT_TRAILING_ZEROS(Int_Time.Next); 
}

/* 
 * Find smallest allowed timestep for local particles 
 * and shared local max & min to these bins
 */

static void set_particle_timebins()
{
	#pragma omp single
	time_bin_min = N_INT_BINS-1, time_bin_max = 0;
	
	#pragma omp for reduction(min:time_bin_min) reduction(max:time_bin_max)
	for (int i = 0; i < NActive_Particles; i++) {
		
		int ipart = Active_Particle_List[i];
		
		float dt = FLT_MAX;

#ifdef GRAVITY
		float dt_cosmo = cosmological_timestep(ipart); 
		
		dt = fmin(dt, dt_cosmo);
#endif
		// add yours here

		dt = fmin(dt, Time.Step_Max);
		
		Assert(dt >= Time.Step_Min, "Timestep too small or not finite ! \n"
				"        ipart=%d, ID=%d, dt=%g, acc=(%g,%g,%g)", 
				ipart, P[ipart].ID, dt, 
				P[ipart].Acc[0], P[ipart].Acc[1], P[ipart].Acc[2]);

		int want = timestep2timebin(dt);

		int allowed = MAX(Time.Max_Active_Bin, P[ipart].Time_Bin);

		P[ipart].Time_Bin = MIN(want, allowed);

		time_bin_min = MIN(time_bin_min, P[ipart].Time_Bin);

		time_bin_max = MAX(time_bin_max, P[ipart].Time_Bin);
	}

	return ;
}

/*
 * set global timestep and handle fullstep. We also have to consider
 * the first and last step separately and stay in sync with the timeline,
 * i.e. we can choose a longer timestep only if it the next time is a 
 * multiple of it.
 */

static void set_global_timestep()
{
	#pragma omp single
	{

	intime_t step_bin = (intime_t) 1 << time_bin_min; // step down

	intime_t step_sync =   
		(intime_t) 1 << COUNT_TRAILING_ZEROS(Int_Time.Current); // stay synced 

	if (Int_Time.Current == Int_Time.Beg) // treat beginning t0
		step_sync = step_bin;
	
	intime_t step_end = Int_Time.End - Int_Time.Current; // don't overstep end

	Int_Time.Step = umin(step_end, umin(step_bin, step_sync));

	Int_Time.Next += Int_Time.Step;

	Time.Next = Integer2Physical_Time(Int_Time.Next);

	Time.Step = Time.Next - Time.Current;

	if (Sig.First_Step) // correct first step
	//if (Int_Time.Current == Int_Time.Beg) // correct first step
		Time.Max_Active_Bin = time_bin_min;

	} // omp single

	#pragma omp flush

	Sig.Fullstep = false;

	if (Int_Time.Current == Int_Time.Next_Full_Step) {
		
		Sig.Fullstep = true;

		#pragma omp single
		Int_Time.Next_Full_Step = Umin(Int_Time.End,  
						Int_Time.Current + ((intime_t) 1 << time_bin_max) );
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

	Time.Step_Max = Time.End - Time.Begin;
	Time.Step_Min =  Time.Step_Max / ((intime_t) 1 << (N_INT_BINS - 1) );

	Time.Max_Active_Bin = N_INT_BINS - 1;

	size_t nBytes = Task.Npart_Total_Max * sizeof(*Active_Particle_List);

	Active_Particle_List = Malloc(nBytes, "Active Part List");

	NActive_Particles = Task.Npart_Total;

	for (int i = 0; i < NActive_Particles; i++) 
		Active_Particle_List[i] = i;

	return ;
}

void Make_Active_Particle_List()
{
	int i = 0;

	#pragma omp single
	{

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		if (P[ipart].Time_Bin <= Time.Max_Active_Bin)
			Active_Particle_List[i++] = ipart; 
	}

	NActive_Particles = i;

	Assert(NActive_Particles > 0, 
			"No Active Particles, instead %d, bin max %d"
			, i, Time.Max_Active_Bin);	
	
	} // omp single

	return ;
}

/* 
 * Give the physical timestep from timebin and vice versa
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
	#pragma omp single nowait
	{

	int npart[N_INT_BINS] = { 0 };
	
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		npart[P[ipart].Time_Bin]++;

	int npart_global[N_INT_BINS] = { 0 };

	MPI_Reduce(npart, npart_global, N_INT_BINS, MPI_INT, MPI_SUM, Sim.Master, 
			MPI_COMM_WORLD);

	if (!Task.Is_MPI_Master)
		goto skip;
	
	int imin = -1, imax = -1;

	for (int i = 0; i < N_INT_BINS; i++)
		if (npart_global[i] != 0 && imin < 0) 
			imin = i;

	for (int i = N_INT_BINS-1; i > -1; i--)
		if (npart_global[i] != 0 && imax < 0)
			imax = i;	

	char fullstep[CHARBUFSIZE] = {" "};

	if (Sig.Fullstep)
		sprintf(fullstep,", Fullstep");

	printf("Systemstep %g, NActive %d %s\n"
			"   Bin       nGas        nDM A    dt\n", 
			Time.Step, NActive_Particles, fullstep );

	for (int i = imax; i > Time.Max_Active_Bin; i--)
		printf("   %2d    %7d     %7d %s  %16.12f \n", 
			i, 0, npart_global[i], " ", Timebin2Timestep(i));

	for (int i = Time.Max_Active_Bin; i >= imin; i--)
		printf("   %2d    %7d     %7d %s  %16.12f \n", 
			i, 0, npart_global[i], "X", Timebin2Timestep(i));

	printf("\n");
	
	if (Sig.Fullstep)
		rprintf("Next full step at t = %g \n\n", 
				Integer2Physical_Time(Int_Time.Next_Full_Step));

	skip:;

	} // omp single nowait

	return ;
}

#undef N_INT_BINS
#undef COUNT_TRAILING_ZEROS

/* 
 * Cosmological N-body step, Dehnen & Read 2011, eq 21
 */

static float cosmological_timestep(const int ipart)
{
	const float acc = ALENGTH3(P[ipart].Acc);
	
	float dt = TIME_INT_ACCURACY * sqrt(2*GRAV_SOFTENING / acc); 
		
	return dt;
}
