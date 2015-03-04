#include "../globals.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

#include "gravity.h"
#include "../domain.h"

#define N_EWALD 64 // Hernquist+ 1991

static inline Float nearest(const Float dx);
static void compute_ewald_force(const int i, const int j, const int k,
								const double r[3], double force[3]);
static double compute_ewald_potential(const double r[3]);

const static double alpha = 2.0;

static double Box2Ewald_Grid = 0;
static double Boxsize = 0;
static double Boxhalf = 0;

static Float Fx[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fy[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fz[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fp[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };

void Gravity_Tree_Periodic()
{

	return ;
}

void Gravity_Tree_Periodic_Init()
{
	rprintf("Init Ewald correction ... ");

	Boxsize = Sim.Boxsize[0];
	Boxhalf = Boxsize / 2;
	Box2Ewald_Grid = 2 * N_EWALD / Boxsize;

	/* this follows closely Gadget-2 */

	int size = p3(N_EWALD + 1) / Sim.NRank; // generic size
	int beg = Task.Rank * size;

	int len = size; // specific size treating last one

	if (Task.Rank == (Sim.NRank - 1))
		len = p3(N_EWALD + 1) - beg;

	int end = beg + len;

	#pragma omp parallel for
	for (int i = 0; i < N_EWALD+1; i++) {

		for (int j = 0; j < N_EWALD+1; j++) {

			for (int k = 0; k < N_EWALD+1; k++) {

				int n = (i * (N_EWALD+1) + j) * (N_EWALD+1) + k;

				if (n < beg || n > end) // on other ranks
					continue;

				double r[3] = { (0.5 * i) / N_EWALD,
								(0.5 * j) / N_EWALD,
								(0.5 * k) / N_EWALD};

				double force[3] = { 0 };

				compute_ewald_force(i, j, k, r, &force[0]);

				Fx[i][j][k] = force[0];
				Fy[i][j][k] = force[1];
				Fz[i][j][k] = force[2];

				Fp[i][j][k] = compute_ewald_potential(r);
			} // k
		} // j
	} // i

	#pragma omp parallel for schedule(static,1)
	for (int rank = 0; rank < Sim.NRank; rank++) {  // merge

		beg = rank * size;

		len = size;

		if (Task.Rank == (Sim.NRank - 1))
			len = p3(N_EWALD + 1) - beg;

		MPI_Bcast(&Fx[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fy[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fz[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fp[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);

	} // for rank

	#pragma omp parallel for
	for (int i = 0; i < N_EWALD+1; i++) {

		for (int j = 0; j < N_EWALD+1; j++) {

			for (int k = 0; k < N_EWALD+1; k++) {

				Fx[i][j][k] /= p2(Boxsize);
				Fy[i][j][k] /= p2(Boxsize);
				Fz[i][j][k] /= p2(Boxsize);
				Fp[i][j][k] /= Boxsize;
			}
		}
	}

	rprintf("done \n");

	return ;
}

static void compute_ewald_correction_table()
{
	return ;
}

const static char fname[CHARBUFSIZE] = { "ewald_table.dat" };

static void read_ewald_correction_table(FILE * fp)
{

	return ;
}

static void write_ewald_corection_table()
{
	const char fname[CHARBUFSIZE] = { "ewald_table.dat" };



	return ;
}

/*
 * Get Ewald correction from grid using CIC binning (Hockney & Eastwood)
 */

static void ewald_correction(const double x, const double y, const double z,
							 double f[3])
{
	f[0] = f[1] = f[2] = 0;

	double qx = x;
	int signx = -1;

	if (qx < 0) {
	
		qx = -qx;
		signx = 1;
	}

	double qy = y;
	int signy = -1;

	if (qy < 0) {
	
		qy = -qy;
		signy = 1;
	}

	double qz = z;
	int signz = -1;

	if (qz < 0) {
	
		qz = -qz;
		signz = 1;
	}

	double u = qx * Box2Ewald_Grid;
	double v = qy * Box2Ewald_Grid;
	double w = qz * Box2Ewald_Grid;

	int i = (int) u;
	int j = (int) v;
	int k = (int) w;

	double weights[8] = { (1 - u) * (1 - v) * (1 - w), // CIC
						  (1 - u) * (1 - v) * w,
						  (1 - u) * v * (1 - w),
						  (1 - u) * v * (w),
						  u * (1 - v) * (1 - w),
						  u * (1 - v) * w,
						  u * v * (1 - w),
						  u * v * w };
	
	f[0] = weights[0] * Fx[i][j][k] + weights[1] * Fx[i][j][k+1] + 
		   weights[2] * Fx[i][j+1][k] + weights[3] * Fx[i][j+1][k+1] + 
		   weights[4] * Fx[i+1][j][k] + weights[5] * Fx[i+1][j][k+1] +
		   weights[6] * Fx[i+1][j+1][k] + weights[7] * Fx[i+1][j+1][k+1];
	f[0] *= signx;

	f[1] = weights[0] * Fy[i][j][k] + weights[1] * Fy[i][j][k+1] + 
		   weights[2] * Fy[i][j+1][k] + weights[3] * Fy[i][j+1][k+1] + 
		   weights[4] * Fy[i+1][j][k] + weights[5] * Fy[i+1][j][k+1] +
		   weights[6] * Fy[i+1][j+1][k] + weights[7] * Fy[i+1][j+1][k+1];
	f[1] *= signy;

	f[2] = weights[0] * Fz[i][j][k] + weights[1] * Fz[i][j][k+1] + 
		   weights[2] * Fz[i][j+1][k] + weights[3] * Fz[i][j+1][k+1] + 
		   weights[4] * Fz[i+1][j][k] + weights[5] * Fz[i+1][j][k+1] +
		   weights[6] * Fz[i+1][j+1][k] + weights[7] * Fz[i+1][j+1][k+1];
	f[2] *= signz;

	return ;
}	

static void ewald_potential(const double x, const double y, const double z,
							double p[1])
{
	p[0] = 0;

	double qx = x;

	if (qx < 0) 
		qx = -qx;

	double qy = y;

	if (qy < 0) 
		qy = -qy;

	double qz = z;

	if (qz < 0) 
		qz = -qz;

	double u = qx * Box2Ewald_Grid;
	double v = qy * Box2Ewald_Grid;
	double w = qz * Box2Ewald_Grid;

	int i = (int) u;
	int j = (int) v;
	int k = (int) w;

	double weights[8] = { (1 - u) * (1 - v) * (1 - w), // CIC
						  (1 - u) * (1 - v) * w,
						  (1 - u) * v * (1 - w),
						  (1 - u) * v * (w),
						  u * (1 - v) * (1 - w),
						  u * (1 - v) * w,
						  u * v * (1 - w),
						  u * v * w };
	
	p[0] = weights[0] * Fp[i][j][k] + weights[1] * Fp[i][j][k+1] + 
		   weights[2] * Fp[i][j+1][k] + weights[3] * Fp[i][j+1][k+1] + 
		   weights[4] * Fp[i+1][j][k] + weights[5] * Fp[i+1][j][k+1] +
		   weights[6] * Fp[i+1][j+1][k] + weights[7] * Fp[i+1][j+1][k+1];

	return ;
}	

/*
 * Compute force correction term from Hernquist+ 1992 (2.14b, 2.16), this
 * is the difference between the force from the infinite lattice and the
 * nearest image. "r" and "force" are in boxsize units.
 */

static void compute_ewald_force(const int i, const int j, const int k,
								const double x[3], double force[3])
{
	force[0] = force[1] = force[2] = 0;
  	
	if(i == 0 && j == 0 && k == 0)
    	return;

	double r = ALENGTH3(x);

	force[0] = x[0] / p3(r);
	force[1] = x[1] / p3(r);
	force[2] = x[2] / p3(r);

	int n[3] = { 0 };

	for (n[0] = -4; n[0] < 5; n[0]++) {

		for (n[1] = -4; n[1] < 5; n[1]++) {
		
			for (n[2] = -4; n[2] < 5; n[2]++) {

				double dx[3] = { x[0]-n[0], x[1]-n[1], x[2]-n[2] };

				double r = ALENGTH3(dx);

				double val = erfc(alpha * r) + 2*alpha * r/sqrt(PI)
							* exp(-p2(alpha*r));

				force[0] -= dx[0] / p3(r) * val;
				force[1] -= dx[1] / p3(r) * val;
				force[2] -= dx[2] / p3(r) * val;
			}
		}
	}

	int h[3] = { 0 };

	for (h[0] = -4; h[0] < 5; h[0]++) {
	
		for (h[1] = -4; h[1] < 5; h[1]++) {
		
			for (h[2] = -4; h[2] < 5; h[2]++) {

				double hdotx = x[0]*h[0] + x[1]*h[1] + x[2]*h[2];

				int h2 = ASCALPROD3(h);

				if (h2 <= 0)
					continue;

				double val = 2.0/h2 * exp(-PI*PI*h2/p2(alpha))
									* sin(2*PI*hdotx);

				force[0] -= h[0] * val;
				force[1] -= h[1] * val;
				force[2] -= h[2] * val;
			}
		}
	}

	return ;
}

static double compute_ewald_potential(const double r[3])
{
	if ((r[0] == 0) && (r[1] == 0) && (r[2] == 0))
		return 2.8372975; // == U, eq. 2.15

	double sum1 = 0;
	int n[3] = { 0 };

	for (n[0] = -4; n[0] < 5; n[0]++) {

		for (n[1] = -4; n[1] < 5; n[1]++) {
		
			for (n[2] = -4; n[2] < 5; n[2]++) {

				double dx[3] = { r[0]-n[0], r[1]-n[1], r[2]-n[2] };

				double dr = ALENGTH3(dx);

				sum1 += erfc(alpha * dr)/dr;
			}
		}
	}

	double sum2 = 0;
	int h[3] = { 0 };

	for (h[0] = -4; h[0] < 5; h[0]++) {

		for (h[1] = -4; h[1] < 5; h[1]++) {
		
			for (h[2] = -4; h[2] < 5; h[2]++) {

				double hdotr = r[0]*h[0] + r[1]*h[1] + r[2]*h[2];

				int h2 = ASCALPROD3(h);

				if (h2 <= 0)
					continue;

				sum2 += 1.0/(PI*h2) * exp(-PI*PI*h2 / p2(alpha)) 
					* cos(2*PI*hdotr);
			}
		}
	}

	return PI/p2(alpha) - sum1 - sum2 + 1.0/ALENGTH3(r);
}

static inline Float nearest(const Float dx)
{
	Float dx_periodic = dx;

	if (dx > Boxhalf)
		dx_periodic = dx - Boxsize;

	if (dx < -Boxhalf)
		dx_periodic = dx + Boxsize;

	return dx_periodic;
}


#undef N_EWALD

#endif
