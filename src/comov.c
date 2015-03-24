#include "globals.h"
#include "timestep.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#ifdef COMOVING

#define TABLESIZE 128

static double Drift_Table[TABLESIZE] = { 0 }, Kick_Table[TABLESIZE] = { 0 },
			  Exp_Factor_Table[TABLESIZE] = { 0 };

static gsl_spline *Kick_Spline = NULL, *Drift_Spline = NULL;
static gsl_interp_accel *Acc[4] = { NULL };
#pragma omp threadprivate(Kick_Spline,Drift_Spline,Acc)

/*
 * These functions get the kick & drift factors of the symplectic integrator
 * in comoving coordinates from a spline interpolation of the integral in
 * Appendix of Quinn+ 1997. This allows us to use a short table resulting in
 * relatives errors < 1e-4, which is the same as the numerical integrator 
 * gives. These functions are thread safe.
 */

double Particle_Kick_Step(const int ipart, const double a_next)
{
	double a_curr = Integer2Physical_Time(P[ipart].Int_Time_Pos);

	double a_beg = gsl_spline_eval(Kick_Spline, a_curr, Acc[0]);
	double a_end = gsl_spline_eval(Kick_Spline, a_next, Acc[1]);

	return a_end - a_beg;
}

double Particle_Drift_Step(const int ipart, const double a_next)
{
	double a_curr = Integer2Physical_Time(P[ipart].Int_Time_Pos);

	double a_beg = gsl_spline_eval(Drift_Spline, a_curr, Acc[2]);
	double a_end = gsl_spline_eval(Drift_Spline, a_next, Acc[3]);

	return a_end - a_beg;
}

/*
 * Integrate in s = \int dt/a^{-2} so we are conserving canonical momentum.
 * See Quinn, Katz, Stadel & Lake 1997, Peebles 1980, Bertschinger 1999. 
 * Note that we are integrating in "da" so the integrals from Quinns paper 
 * have to be transformed from dt, which gives the additional factor of 
 * 1/\dot{a} = 1/H(a)/a. 
 
 */

static double comoving_symplectic_drift_integrant(double a, void *param)
{
	return 1 / (Hubble_Parameter(a) * a*a*a);
}

static double comoving_symplectic_kick_integrant(double a, void *param)
{
	return 1 / (Hubble_Parameter(a) * a*a);
}

/*
 * This sets up the comoving kick & drift factors using a table and cspline
 * interpolation. For an EdS universe the solution is:
 * drift(a) = -2/H0/sqrt(a) ; kick(a) = 2/H0 * sqrt(a)
 */

void Setup_Comoving()
{
	#pragma omp parallel
	{

	gsl_function gsl_F = { 0 };
	gsl_integration_workspace *gsl_workspace = NULL;
	gsl_workspace = gsl_integration_workspace_alloc(TABLESIZE);

	const double time_min = Time.Begin;
	const double da = log(Time.End) - log(Time.Begin);

	#pragma omp for
	for (int i = 0; i < TABLESIZE; i++) {

		double error = 0;

		double time_max = exp(log(Time.Begin) + da * i/(TABLESIZE - 1.0) );

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

	} // omp parallel

	return ;
}

void Finish_Comoving()
{
	Free(Kick_Table);
	Free(Drift_Table);
	Free(Exp_Factor_Table);

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


#endif // COMOVING
