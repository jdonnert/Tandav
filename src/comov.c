#include "globals.h"
#include "timestep.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#ifdef COMOVING

#define TABLESIZE 128

#ifdef GADGET_COMOVING_VEL_UNIT
static void convert_velocities_to_comoving();
#else
static inline void convert_velocities_to_comoving(){};
#endif


static void setup_kick_drift_factors();

static double Drift_Table[TABLESIZE] = { 0 }, Kick_Table[TABLESIZE] = { 0 },
			  Exp_Factor_Table[TABLESIZE] = { 0 };

static gsl_spline *Kick_Spline = NULL, *Drift_Spline = NULL;
static gsl_interp_accel *Acc[4] = { NULL };
#pragma omp threadprivate(Kick_Spline,Drift_Spline,Acc)

/*
 * These functions get the kick & drift factors of the symplectic integrator
 * in comoving coordinates from a spline interpolation of the integral in
 * Appendix of Quinn+ 1997. This allows us to use a short table resulting in
 * relative errors < 1e-4, which is the same as the numerical integrator 
 * gives. As our velocity variable is \dot{x}*a^1.5, but the canonical momentum
 * is m*a^2\dot{x}, we have to add sqrt(a) to the particle drift step.
 * These functions are thread safe.
 */

double Particle_Kick_Step(const int ipart, const double a_next)
{
	double a_curr = Integer_Time2Integration_Time(P.Int_Time_Pos[ipart]);

	double kick_factor_beg = gsl_spline_eval(Kick_Spline, a_curr, Acc[0]);
	double kick_factor_end = gsl_spline_eval(Kick_Spline, a_next, Acc[1]);

	return kick_factor_end - kick_factor_beg;
}

double Particle_Drift_Step(const int ipart, const double a_next)
{
	double a_curr = Integer_Time2Integration_Time(P.Int_Time_Pos[ipart]);

	double drift_factor_beg = gsl_spline_eval(Drift_Spline, a_curr, Acc[2]);
	double drift_factor_end = gsl_spline_eval(Drift_Spline, a_next, Acc[3]);

	return (drift_factor_end - drift_factor_beg);
}

/*
 * Integrate in s = \int^{t_i}_{t_0} dt/a^{-2} so we are conserving 
 * canonical momentum m * a^2 * \dot{x}.
 * See Quinn, Katz, Stadel & Lake 1997, Peebles 1980, Bertschinger 1999. 
 * Note that in code units "dt = da" so the integrals from Quinns paper 
 * have to be transformed from dt -> da, which gives the additional factor of 
 * 1/\dot{a} = 1/H(a)/a. Time stepping however is done in dln(a) = 1+z.
 */

static double comoving_symplectic_drift_integrant(double a, void *param)
{
	return 1.0 / (Hubble_Parameter(a) * a*a*a);
}

static double comoving_symplectic_kick_integrant(double a, void *param)
{
	return 1.0 / (Hubble_Parameter(a) * a*a);
}

void Setup_Comoving()
{
	#pragma omp parallel
	{

	convert_velocities_to_comoving(); // GADGET_COMOVING_VEL_UNIT

    setup_kick_drift_factors();

	} // omp parallel

	return ;
}

void Finish_Comoving()
{
	#pragma omp parallel
	{

	gsl_spline_free(Drift_Spline);
	gsl_spline_free(Kick_Spline);

	gsl_interp_accel_free(Acc[0]);
	gsl_interp_accel_free(Acc[1]);
	gsl_interp_accel_free(Acc[2]);
	gsl_interp_accel_free(Acc[3]);

	} // omp parallel

	return ;
}

/*
 * This sets up the comoving kick & drift factors using a table and cspline
 * interpolation. For an EdS universe the solution is:
 * drift(a) = -2/H0/sqrt(a) ; kick(a) = 2/H0 * sqrt(a)
 */

static void setup_kick_drift_factors()
{
	gsl_function gsl_F = { 0 };
	gsl_integration_workspace *gsl_workspace = NULL;
	gsl_workspace = gsl_integration_workspace_alloc(TABLESIZE);

	const double time_min = Time.Begin;
	const double da = (log(Time.End) - log(Time.Begin)) /(TABLESIZE - 1.0);

	#pragma omp for
	for (int i = 0; i < TABLESIZE; i++) {

		double error = 0;

		double time_max = exp(log(Time.Begin) + da * i );

		Exp_Factor_Table[i] = time_max;

		gsl_F.function = &comoving_symplectic_drift_integrant;

		gsl_integration_qag(&gsl_F, time_min, time_max, 0, 1e-8, TABLESIZE,
				GSL_INTEG_GAUSS41, gsl_workspace, &Drift_Table[i], &error);

		gsl_F.function = &comoving_symplectic_kick_integrant;

		gsl_integration_qag(&gsl_F, time_min, time_max, 0, 1e-8, TABLESIZE,
				GSL_INTEG_GAUSS41, gsl_workspace, &Kick_Table[i], &error);

	} // for i

	gsl_integration_workspace_free(gsl_workspace);

	Acc[0] = gsl_interp_accel_alloc();
	Acc[1] = gsl_interp_accel_alloc();
	Acc[2] = gsl_interp_accel_alloc();
	Acc[3] = gsl_interp_accel_alloc();

	Drift_Spline = gsl_spline_alloc(gsl_interp_cspline, TABLESIZE);
	Kick_Spline = gsl_spline_alloc(gsl_interp_cspline, TABLESIZE);

	gsl_spline_init(Kick_Spline, Exp_Factor_Table, Kick_Table, TABLESIZE);
	gsl_spline_init(Drift_Spline, Exp_Factor_Table, Drift_Table, TABLESIZE);

	return;
}

/*
 * This converts the particle velocities from the initial conditions to the 
 * internal velocity variable u = v*a^1.5. This is because in an Einstein-
 * de Sitter universe the peculiar velocity at high redshift scales with
 * a^0.5, as can be shown by combining the Newtonian equations of motions in
 * comoving coordinates with the Zeldovich approximation. However, we use a 
 * modified comoving potential from Quinn+ 1997 to construct the symplectic 
 * integrator, which requires another factor "a" in the comoving potential 
 * and hence the peculiar velocity scales with a^1.5 using this comoving 
 * potential. (Peebles 1980, Mo, v.d.Bosch, White 4.1.8, Springel 2001)
 */

#ifdef GADGET_COMOVING_VEL_UNIT
static void convert_velocities_to_comoving()
{
	const double phys2comov_vel = pow(Time.Current, 1.5);

	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		P.Vel[0][ipart] *= phys2comov_vel;
		P.Vel[1][ipart] *= phys2comov_vel;
		P.Vel[2][ipart] *= phys2comov_vel;
	}

	return ;
}
#endif

/* 
 * In cosmological simulations the timestep has to be bound by the maximum 
 * displacement given the rms velocity of all particles. The maximum 
 * displacement is set relative to the mean particle separation.
 * (Knebe et al. 2009, Lukic et al. 2007, Crocce et al. 2006)
 */

static double vel2[NPARTYPE] = { 0 }, min_mpart[NPARTYPE] = { 0 };
static long long npart[NPARTYPE] = { 0 };

double Comoving_VelDisp_Timestep_Constraint(const double dt_max_ext)
{
	#pragma omp single
	{

	for (int type = 0; type < NPARTYPE; type++) {

		vel2[type] = npart[type] = 0;
		min_mpart[type] = DBL_MAX;
	}

	} // omp single

	double vel2_thread[NPARTYPE] = { 0 };
	double min_mpart_thread[NPARTYPE] = { DBL_MAX };
	int npart_thread[NPARTYPE] = { 0 };

	#pragma omp for nowait
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		int type = P.Type[ipart];

		vel2_thread[type] += p2(P.Vel[0][ipart]) + p2(P.Vel[1][ipart]) 
							+ p2(P.Vel[2][ipart]);

		min_mpart_thread[type] = fmin(min_mpart[type], P.Mass[ipart]);

		npart_thread[type]++;
	}

	#pragma omp critical // array-reduce 
	{

	for (int type = 0; type < NPARTYPE; type++) {

		vel2[type] += vel2_thread[type];

		min_mpart[type] = fmin(min_mpart_thread[type], min_mpart[type]);

		npart[type] += npart_thread[type];
	} // for type

	} // omp critical

	#pragma omp barrier

	double dt_max = DBL_MAX;

	#pragma omp single copyprivate(dt_max)
	{

	MPI_Allreduce(MPI_IN_PLACE, vel2, NPARTYPE, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, min_mpart, NPARTYPE, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, npart, NPARTYPE, MPI_LONG_LONG, MPI_SUM,
			MPI_COMM_WORLD);

	double rho_baryon = Cosmo.Omega_Baryon*Cosmo.Rho_Crit0;
	double rho_nonbaryon = Cosmo.Omega_Matter*Cosmo.Rho_Crit0 - rho_baryon;

	printf("\nComoving Time Displacement Constraint at a = %g \n",
		   Time.Current);

	for (int type = 0; type < NPARTYPE; type++) {

		if (npart[type] == 0)
			continue;

		double dmean = 0;

		if (type == 0)
			dmean = pow( min_mpart[type]/rho_baryon, 1.0/3.0);
		else
			dmean = pow( min_mpart[type]/rho_nonbaryon, 1.0/3.0);

		double vrms = sqrt( vel2[type]/npart[type] );

		double dt = TIME_DISPL_CONSTRAINT * dmean / vrms
					* p2(Cosmo.Expansion_Factor) * Cosmo.Hubble_Parameter ;

		printf("   Type %d: dmean %g, min mpart %g, sqrt(<p^2>) %g,"
				" dlogmax %g \n", type, dmean, min_mpart[type], vrms, dt);

		dt_max = fmin(dt, dt_max);
	}

	} // omp single

	return fmin(dt_max_ext, dt_max);
}

#endif // COMOVING
