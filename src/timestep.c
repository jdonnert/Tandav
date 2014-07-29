#include "globals.h"
#include "timestep.h"

static int timestep2timebin(const float dt);
static void print_timebins();

static float cosmological_timestep(const int ipart);
static float kepler_timestep(const int ipart);
static float pseudo_symmetric_kepler_timestep(const int ipart);

struct TimeData Time;


void Set_New_Timesteps()
{
	Profile("Timesteps");

	int local_bin_min = 65, local_bin_max = 0;

	#pragma omp parallel for \
		reduction(min:local_bin_min) reduction(max:local_bin_max)
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {

		double dt = DBL_MAX;

		double dt_cosmo = cosmological_timestep(ipart); 
		
		//dt = fmin(dt, dt_cosmo);
		
		//dt = 8.56948 / 1000.;
if(ipart == 0)
	dt=8;
		double dt_kepler = kepler_timestep(ipart);
		
		dt = fmin(dt, dt_kepler);
		
		/* Add above */
		Assert(dt > Time.StepMin, 
				"Timestep too small ! Ipart=%d, dt=%g, dt_cosmo=%g",
				ipart, dt, dt_cosmo);

		P[ipart].TimeBin = timestep2timebin(dt);
		
		local_bin_min = Imin(local_bin_min, P[ipart].TimeBin);
		local_bin_max = Imax(local_bin_max, P[ipart].TimeBin);
	}
	
	int global_bin_min = 65;
	
	MPI_Allreduce(&local_bin_min, &global_bin_min, 1, MPI_INT, MPI_MIN, 
			MPI_COMM_WORLD);

	Time.IntStep = 1ULL << global_bin_min;
	Time.Step = Timebin2Timestep(global_bin_min);
	
	Flag.Fullstep = false;

	print_timebins();

	if (Time.IntCurrent == Time.IntSyncPoint) {
		
		Flag.Fullstep = true;

		int global_bin_max = 0;

		MPI_Allreduce(&local_bin_max, &global_bin_max, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);
		
		Time.IntSyncPoint = Time.IntCurrent + 1ULL << global_bin_max;
	
		rprintf("Next sync point at %g \n", 
				Integer2PhysicalTime(Time.IntSyncPoint));
	}

	Profile("Timesteps");
Flag.Endrun = true;
	return ;
}

/* Cosmological N-body step, Dehnen & Read 2011, eq (21) */
static float cosmological_timestep(const int ipart)
{
	const float acc = len3(P[ipart].Force) / P[ipart].Mass;
	const float surge = P[ipart].Surge;

	float dt = sqrt(2*TIME_INT_ACCURACY * GRAV_SOFTENING / acc); 
		
//printf("dt=%g surge=%g acc=%g term=%g \n", dt, surge, acc, - dt  * surge / acc);
	float deriv = 0.5 * dt  * surge  /acc;

	if (Time.Current != Time.Begin)
		dt /= (1 + deriv );

	return dt;
}

/* Timestep for the Kepler problem,  Dehnen & Read 2011, (20)*/
static float kepler_timestep(const int ipart)
{
	double dx = P[0].Pos[0] - P[1].Pos[0];
	double dy = P[0].Pos[1] - P[1].Pos[1];
	double dz = P[0].Pos[2] - P[1].Pos[2];

	double r = sqrt( dx*dx + dy*dy + dz*dz );

	float phi = Const.Gravity * P[ipart].Mass / r;

 	float acc_phys = len3(P[ipart].Force) / P[ipart].Mass;

	float dt_kepler = TIME_INT_ACCURACY * sqrt(phi) / acc_phys;

	return dt_kepler;
}

static float pseudo_symmetric_kepler_timestep(const int ipart)
{
	double dx = P[0].Pos[0] - P[1].Pos[0];
	double dy = P[0].Pos[1] - P[1].Pos[1];
	double dz = P[0].Pos[2] - P[1].Pos[2];

	double r = sqrt( dx*dx + dy*dy + dz*dz );

	double mu = P[0].Mass * Const.Gravity;

 	float acc_phys = len3(P[ipart].Force) / P[ipart].Mass;

	float dt = TIME_INT_ACCURACY *  sqrt(r)/acc_phys;

	double dvx = P[0].Vel[0] - P[1].Vel[0];
	double dvy = P[0].Vel[1] - P[1].Vel[1];
	double dvz = P[0].Vel[2] - P[1].Vel[2];

	float vdotr = sqrt( dvx*dx + dvy*dy + dvz*dz );

	float deriv = 1.5 * dt * vdotr / r/r; 

	if (Time.Current != Time.Begin)
		dt /= (1 - 0.5 * deriv );

	return dt;
}


void Setup_Timesteps()
{
	Time.IntBeg = 0;
	Time.IntEnd = 0xFFFFFFFFFFFFFFFF;
	Time.IntCurrent = Time.IntBeg;

	Time.StepMax = Time.End - Time.Begin;
	Time.StepMin = Time.StepMax / (1ULL << 63);

	return ;
}

bool Time_Is_Up()
{
	if (Time.IntCurrent == Time.IntEnd)
		return true;

	if (Flag.Endrun)
		return true;

	if (Runtime() >= Param.RuntimeLimit) {

		Flag.WriteRestartFile = true;

		return true;
	}

	return false;
}

bool Time_For_Snapshot()
{
	if (Flag.WriteSnapshot)
		return true;
 	
	if (Time.Current + Time.Step > Time.NextSnap) {
	
		Time.NextSnap += Time.BetSnap;

		rprintf("Snapshot at t=%g, Next Snapshot at t=%g \n", 
				Time.Current, Time.NextSnap);

		return true;
	}

	return false;
}

/* Give the physical timestep from timebin
 * Use Time.Step for the current system step */
double Timebin2Timestep(const int TimeBin)
{
	return Time.StepMin * (1ULL << TimeBin);
}


double Integer2PhysicalTime(const uint64_t IntTime)
{
	return Time.Begin +  IntTime * Time.StepMin;
}

/* Convert a timestep to a power 2 based timebin 
 * it is assumed that dt >= StepMin */
static int timestep2timebin(const double dt)
{
	return 63 - ceil( log2( (Time.End - Time.Begin)/dt ) );
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

	rprintf("Timesteps: shortest = %d , longest = %d\n"
			"   Bin       nGas           nDM         dt\n", imax, imin);

	for (int i = imax; i >= imin; i--)
			rprintf("   %2d    %7d        %7d       %g \n", 
					i, 0, npart_global[i], Timebin2Timestep(i));


	rprintf("\n");

	MPI_Barrier(MPI_COMM_WORLD);

	return ;
}

