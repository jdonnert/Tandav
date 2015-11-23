#include "globals.h"
#include "timestep.h"

/* 
 * The number of bins is given by the number of bits in an integer time 
 */

#define N_INT_BINS (sizeof(intime_t) * CHAR_BIT)
#define COUNT_TRAILING_ZEROS(x) __builtin_ctzll(x)

static int max_active_time_bin();
static void set_new_particle_timebins();
static void set_system_timestep();
static int timestep2timebin(const double dt);
static void print_timebins();

static void set_global_timestep_constraint();
static float get_physical_timestep(const int);
static float cosmological_timestep(const int ipart, const Float accel);

struct TimeData Time = { 0 };
struct IntegerTimeLine Int_Time = { 0 };

static float Dt_Max_Global = FLT_MAX;
static int Time_Bin_Min = N_INT_BINS-1, Time_Bin_Max = 0;

/* 
 * All active particles get a new step that is smaller than or equal to 
 * the largest active bin. We also set the fullstep signal. 
 */

void Set_New_Timesteps()
{
	Profile("Timesteps");

	if (Sig.Sync_Point || Sig.First_Step)
		set_global_timestep_constraint(); // set t_max_global

	set_new_particle_timebins();

	#pragma omp single
	{

	MPI_Allreduce(MPI_IN_PLACE, &Time_Bin_Min, 1, MPI_INT, MPI_MIN,
			MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, &Time_Bin_Max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

	set_system_timestep(Time_Bin_Max, Time_Bin_Min);

	Time.Max_Active_Bin = max_active_time_bin();

	} // omp single

	Sig.Sync_Point = false;

	if ((Int_Time.Current == Int_Time.Next_Sync_Point)) {

		Sig.Sync_Point = true;

		#pragma omp barrier

		#pragma omp single
		Int_Time.Next_Sync_Point += 1ULL << Time_Bin_Max;
	}

	Make_Active_Particle_List();

	#pragma omp master
	{

	print_timebins();

	Warn(Time.Step < Param.Min_Timestep, "Time step %g has fallen below "
			"Min_Timestep parameter %g", Time.Step, Param.Min_Timestep);

	} // omp master

	#pragma omp barrier

	Profile("Timesteps");

	return ;
}

/* 
 * The timeline is represented by an integer, where an increment of one 
 * corresponds to the whole integration time divided by 2^(N_INT_BINS-1).
 * In comoving coordinates the timesteps are divided in log space, which
 * means, we are effectively stepping in redshift. Time integration is in da.
 */

void Setup_Time_Integration()
{
	Time.Next_Snap = Time.First_Snap;

	Time.NSnap = (Time.End - Time.Begin)/Time.Bet_Snap + 1;

	rprintf("\nSimulation timeline: \n"
			"   start = %g, end = %g, delta = %g, NSnap = %d \n",
			Time.Begin, Time.End, Time.Bet_Snap, Time.NSnap);

#ifdef COMOVING
	rprintf("   initial redshift = %g, final redshift = %g \n",
			1.0/Time.Begin - 1, 1.0/Time.End - 1);
#endif // COMOVING

	rprintf("\n");

	Assert(Time.NSnap > 0, "Timeline does not seem to produce any outputs");

	Int_Time.Beg = (intime_t) 0;
	Int_Time.End = (intime_t) 1 << (N_INT_BINS - 1);
	Int_Time.Current = Int_Time.Beg;

#ifdef COMOVING
	Time.Step_Max = log(Time.End) - log(Time.Begin); // step in log(a)
#else
	Time.Step_Max = Time.End - Time.Begin;
#endif // ! COMOVING

	Time.Step_Min = Time.Step_Max / Int_Time.End;

	Time.Current = Integer_Time2Integration_Time(Int_Time.Beg);

	Time.Max_Active_Bin = N_INT_BINS - 1;

	size_t nBytes = Task.Npart_Total_Max * sizeof(*Active_Particle_List);

	Active_Particle_List = Malloc(nBytes, "Active Part List");

	NActive_Particles = Task.Npart_Total;

	#pragma omp parallel for	
	for (int i = 0; i < NActive_Particles; i++) 
		Active_Particle_List[i] = i;

	return ;
}

static inline Float convert_dt_to_dlna(const Float dt)
{
#ifdef COMOVING
	return dt * Cosmo.Hubble_Parameter;
#else 
	return dt;
#endif
}

/* 
 * Find smallest allowed timestep for local particles given the time step 
 * criteria. Find local max & min to these bins. 
 */


static void set_new_particle_timebins()
{
	#pragma omp single
	{

	Time_Bin_Min = N_INT_BINS-1;
	Time_Bin_Max = 0;

	} // omp single

	#pragma omp for reduction(min:Time_Bin_Min) reduction(max:Time_Bin_Max)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		Float dt = get_physical_timestep(ipart);

		dt = convert_dt_to_dlna(dt); // COMOVING

		dt = fmin(dt, Dt_Max_Global);

		Assert(dt >= Time.Step_Min, "Timestep too small for integer timeline"
				" or not finite ! \n        ipart=%d, ID=%d, dt=%g, "
				"acc=(%g,%g,%g)", ipart, P[ipart].ID, dt,
				P[ipart].Acc[0], P[ipart].Acc[1], P[ipart].Acc[2]);

		int want = timestep2timebin(dt);

		int allowed = MAX(Time.Max_Active_Bin, P[ipart].Time_Bin);
		
		P[ipart].Time_Bin = MIN(want, allowed);

		Time_Bin_Min = MIN(Time_Bin_Min, P[ipart].Time_Bin);
		Time_Bin_Max = MAX(Time_Bin_Max, P[ipart].Time_Bin);
	}


	return ;
}

/*
 * Set global system timestep. We also have to consider
 * the first and last step separately and stay in sync with the timeline,
 * i.e. we can choose a longer timestep only if it the next time is a 
 * multiple of it.
 */

static void set_system_timestep()
{
	intime_t step_bin = (intime_t) 1 << Time_Bin_Min; // step down ?

	intime_t step_sync = 1ULL << COUNT_TRAILING_ZEROS(Int_Time.Current);

	if (Int_Time.Current == Int_Time.Beg) // treat beginning t0
		step_sync = step_bin;

	intime_t step_end = Int_Time.End - Int_Time.Current; // don't overstep end

	Int_Time.Step = umin(step_end, umin(step_bin, step_sync));

	Int_Time.Next =  Int_Time.Current + Int_Time.Step;

	Time.Next = Integer_Time2Integration_Time(Int_Time.Next);

	Time.Step = Time.Next - Time.Current; // if COMOVING in a

	if (Sig.First_Step)
		Time.Max_Active_Bin = Time_Bin_Min; // correct first step

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

void Make_Active_Particle_List()
{
	#pragma omp single
	{

	int i = 0;

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		if (P[ipart].Time_Bin <= Time.Max_Active_Bin)
			Active_Particle_List[i++] = ipart;
	}

	NActive_Particles = i;

	Assert(NActive_Particles > 0, "No Active Particles, instead %d, bin max %d"
			, i, Time.Max_Active_Bin);

	} // omp single

	return ;
}

/* 
 * Give the integration timestep from timebin and convert from integer to 
 * integration time. In comoving coordinates/cosmological simulations we 
 * multi-step in "dlog(a) = 1+z. We return here dln(a) from the timebin. 
 * Note that dt = da / H(a).
 */

#ifdef COMOVING

double Timebin2Timestep(const int TimeBin)
{
	return Time.Step_Min*((intime_t) 1 << TimeBin); // in dlog(a)
}

double Integer_Time2Integration_Time(const intime_t Integer_Time)
{
	return Time.Begin * exp(Integer_Time * Time.Step_Min); // a
}

double Integer2Physical_Time(const intime_t Integer_Time)
{
	return Integer_Time2Integration_Time(Integer_Time) / Cosmo.Hubble_Parameter;
}

#else // ! COMOVING

double Timebin2Timestep(const int TimeBin)
{
	return Time.Step_Min * ((intime_t) 1 << TimeBin);
}

double Integer2Physical_Time(const intime_t Integer_Time)
{
	return Time.Begin + Integer_Time * Time.Step_Min;
}

double Integer_Time2Integration_Time(const intime_t Integer_Time)
{
	return Integer2Physical_Time(Integer_Time);
}

#endif // ! COMOVING


/* 
 * Convert a timestep to a power 2 based timebin via ceil(log2(StepMax/dt)) 
 */

static int timestep2timebin(const double dt)
{
	int exponent;

	frexp(Time.Step_Max/dt, &exponent);

	return N_INT_BINS - 1 - exponent;
}

static void print_timebins()
{
	int npart[N_INT_BINS] = { 0 };

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
		npart[P[ipart].Time_Bin]++;

	MPI_Reduce(MPI_IN_PLACE, npart, N_INT_BINS, MPI_INT, MPI_SUM, Sim.Master,
			MPI_COMM_WORLD);

	if (!Task.Is_MPI_Master)
		goto skip;

	int imin = -1, imax = -1;

	for (int i = 0; i < N_INT_BINS; i++)
		if ((npart[i] != 0) ) {

			imin = i;

			break;
		}

	for (int i = N_INT_BINS-1; i > -1; i--)
		if ((npart[i] != 0)) {

			imax = i;

			break;
		}

	char step[CHARBUFSIZE] = {"Point"};

	if (Sig.Sync_Point)
		sprintf(step,"Sync point");

#ifdef COMOVING
	rprintf("\n%s <%d> a = %g -> %g\n\n"
			"NActive %d, z = %g, da = %g \n"
			"   Bin       nGas        nDM A    dlog(a)\n",
			step, Time.Step_Counter, Time.Current,
			Integer_Time2Integration_Time(Int_Time.Next),
			NActive_Particles, 1/Time.Current-1, Time.Step);
#else
	rprintf("\n%s <%d> t = %g -> %g\n\n"
			"Systemstep %g, NActive %d \n"
			"   Bin       nGas        nDM A    dt\n",
			step, Time.Step_Counter, Time.Current,
			Integer_Time2Integration_Time(Int_Time.Next),
			Time.Step, NActive_Particles );
#endif // ! COMOVING

	for (int i = imax; i > Time.Max_Active_Bin; i--)
		printf("   %2d    %7d     %7d %s  %16.12f \n",
				i, 0, npart[i], " ", Timebin2Timestep(i));

	for (int i = MIN(Time.Max_Active_Bin, imax); i >= imin; i--)
		printf("   %2d    %7d     %7d %s  %16.12f \n",
			i, 0, npart[i], "X", Timebin2Timestep(i));

	printf("\n");

	if (Sig.Sync_Point)
		rprintf("Next sync point at t = %g \n\n",
				Integer_Time2Integration_Time(Int_Time.Next_Sync_Point));

	skip:;

	return ;
}

#undef N_INT_BINS
#undef COUNT_TRAILING_ZEROS

/*
 * Collect minimum of all timesteps 
 */

static float get_physical_timestep(const int ipart)
{
	const Float acc_phys = Acceleration_Physical(ipart);

	float dt = FLT_MAX;

#ifdef GRAVITY
	float dt_cosmo = cosmological_timestep(ipart, acc_phys);

	dt = fmin(dt, dt_cosmo);
#endif

	// add yours here

	return dt;
}

/* 
 * Cosmological N-body step, Dehnen & Read 2011, eq 21
 */

static float cosmological_timestep(const int ipart, const Float acc_phys)
{
#ifdef COMOVING
	return TIME_INT_ACCURACY * sqrt(2 * Cosmo.Expansion_Factor
			* GRAV_SOFTENING / acc_phys);
#else
	return TIME_INT_ACCURACY * sqrt(2 * GRAV_SOFTENING / acc_phys);
#endif // ! COMOVING
}

/*
 * Collect timestep upper bounds from various models
 */

static void set_global_timestep_constraint()
{
	double dt = Param.Max_Timestep;

	dt = Comoving_VelDisp_Timestep_Constraint(dt); // COMOVING

	// add yours here

	#pragma omp single
	Dt_Max_Global = dt;

	rprintf("Found max global timestep  %g \n", Dt_Max_Global);

	return ;
}
