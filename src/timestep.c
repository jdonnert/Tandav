#include "globals.h"
#include "timestep.h"

struct TimeData Time;

static uint64_t timestep2timebin(const float dt);
static void print_timebins();

static float cosmological_timestep(const int ipart);
static float kepler_timestep(const int ipart);
static float pseudo_symmetric_kepler_timestep(const int ipart);

void Set_New_Timesteps()
{
	Profile("Timesteps");

	float dt_min = FLT_MAX; 

#pragma omp parallel for reduction(min:dt_min) // OpenMP 3.1
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {

		float dt = Time.StepMax;

		float dt_cosmo = cosmological_timestep(ipart); 
		
		dt = fmin(dt, dt_cosmo);
		
		//dt = 8.56948 / 1000.;

		//float dt_kepler = pseudo_symmetric_kepler_timestep(ipart);

		//dt = dt_kepler;//fmin(dt, dt_kepler);
		
		/* Add above */
		Assert(dt > Time.StepMin, 
				"Timestep too small ! Ipart=%d, dt=%g, dt_cosmo=%g",
				ipart, dt, dt_cosmo);

		P[ipart].TimeBin = timestep2timebin(dt);
		
		dt_min = fmin(dt_min, dt);
	}

	MPI_Allreduce(&dt_min, &Time.Step, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
	
	//print_timebins();

	Profile("Timesteps");

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
	Time.StepMax = Time.End - Time.Begin;

	Time.StepMin = Time.StepMax / ((uint64_t) 1 << 63);

	return ;
}

bool Time_Is_Up()
{
	if (Time.Current + Time.Step > Time.End)
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

/* Give the physical timestep of particle ipart 
 * Use Time.Step for the system step */
float Timestep(const int ipart)
{
	return Time.StepMin * (1 << P[ipart].TimeBin);
}

/* Convert a timestep to a power 2 based timebin 
 * it is assumed that dt >= StepMin */
static uint64_t timestep2timebin(const float dt)
{
	return floor( log2(dt / Time.StepMin) );
}

static void print_timebins()
{
	int npart[64] = { 0 };
	
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++)
		npart[P[ipart].TimeBin]++;

	int npart_global[64] = { 0 };

	MPI_Reduce(npart, npart_global, 64, MPI_SUM, MPI_INT, MASTER, 
			MPI_COMM_WORLD);

	rprintf("");
	rprintf("");
	
	return ;
}

