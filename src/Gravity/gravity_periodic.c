#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_periodic.h"

#if defined(GRAVITY) && defined (PERIODIC)

#define N_EWALD 64 // Hernquist+ 1991

static void compute_ewald_correction_table();
static void write_ewald_correction_table();
static bool read_ewald_correction_table();
static void compute_ewald_force(const int, const int, const int,
								const double r[3], double force[3]);
static double compute_ewald_potential(const double r[3]);

const static double Alpha = 2.0; // Hernquist+ 1991 (2.13)
const static char fname[CHARBUFSIZE] = { "./ewald_tables.dat" };

static double Box2Ewald_Grid = 0;

static Float Fx[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { { { 0 } } };
static Float Fy[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { { { 0 } } };
static Float Fz[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { { { 0 } } };
static Float Fp[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { { { 0 } } };

/*
 * Get Ewald correction from the grid using a modified CIC binning 
 * (Hockney & Eastwood) to exploit the symmetry of the Ewald correction in 
 * the octants of the grid. This way the grid of size N_EWALD is effectively
 * doubled in resolution.
 */

void Ewald_Correction(const Float dr[3], Float f[3])
{
	int sign[3] = { -1, -1, -1 };
	
	double dx = dr[0];

	if (dx < 0) {

		dx = -dx;
		sign[0] = 1;
	}

	double dy = dr[1];

	if (dy < 0) {

		dy = -dy;
		sign[1] = 1;
	}

	double dz = dr[2];

	if (dz < 0) {

		dz = -dz;
		sign[2] = 1;
	}

	double u = dx * Box2Ewald_Grid;
	double v = dy * Box2Ewald_Grid;
	double w = dz * Box2Ewald_Grid;

	int i = (int) u;
	int j = (int) v;
	int k = (int) w;

	if (i >= N_EWALD) // don't overshoot
		i = N_EWALD-1;

	if (j >= N_EWALD)
		j = N_EWALD-1;

	if (k >= N_EWALD)
		k = N_EWALD-1;

	u -= i;
	v -= j;
	w -= k;

	double weights[8] = { (1 - u) * (1 - v) * (1 - w), // CIC with u,v,w < 2 !
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

	f[1] = weights[0] * Fy[i][j][k] + weights[1] * Fy[i][j][k+1] +
		   weights[2] * Fy[i][j+1][k] + weights[3] * Fy[i][j+1][k+1] +
		   weights[4] * Fy[i+1][j][k] + weights[5] * Fy[i+1][j][k+1] +
		   weights[6] * Fy[i+1][j+1][k] + weights[7] * Fy[i+1][j+1][k+1];

	f[2] = weights[0] * Fz[i][j][k] + weights[1] * Fz[i][j][k+1] +
		   weights[2] * Fz[i][j+1][k] + weights[3] * Fz[i][j+1][k+1] +
		   weights[4] * Fz[i+1][j][k] + weights[5] * Fz[i+1][j][k+1] +
		   weights[6] * Fz[i+1][j+1][k] + weights[7] * Fz[i+1][j+1][k+1];

	f[0] *= sign[0];
	f[1] *= sign[1];
	f[2] *= sign[2];

	return ;
}

#ifdef GRAVITY_POTENTIAL
void Ewald_Potential(const double dr[3], Float p[1])
{
	double dx = dr[0];

	if (dx < 0)
		dx = -dx;

	double dy = dr[1];

	if (dy < 0)
		dy = -dy;

	double dz = dr[2];

	if (dz < 0)
		dz = -dz;

	double u = dx * Box2Ewald_Grid;
	double v = dy * Box2Ewald_Grid;
	double w = dz * Box2Ewald_Grid;

	int i = (int) u;
	int j = (int) v;
	int k = (int) w;

	if (i >= N_EWALD)
		i = N_EWALD-1;

	if (j >= N_EWALD)
		j = N_EWALD-1;

	if (k >= N_EWALD)
		k = N_EWALD-1;

	u -= i;
	v -= j;
	w -= k;

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
#endif // GRAVITY_POTENTIAL

/*
 * This initialises the cubes holding the Ewald correction force and
 * potential following Hernquist, Bouchet & Suto 1991. Most of the 
 * code is shamelessly copied from Gadget-2 (Springel 2006).
 */

void Gravity_Periodic_Init()
{
	rprintf("Init Ewald correction ");

	Warn((Sim.Boxsize[0] != Sim.Boxsize[1])
		&& (Sim.Boxsize[0] != Sim.Boxsize[2]),
		 "Periodic Gravity requires cubic boxes.\n"
		 "                Setting global Boxsize %g",
		 Sim.Boxsize[0]);

	Boxsize = Sim.Boxsize[1] = Sim.Boxsize[2] = Sim.Boxsize[0];

	Boxhalf = Boxsize / 2;

	Box2Ewald_Grid = 2 * N_EWALD / Boxsize; // clever symmetry mapping

	bool table_found = false;

	if (Task.Is_Master)
		 table_found = read_ewald_correction_table();

	MPI_Bcast(&table_found, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	if (table_found) {
	
		int len = p3(N_EWALD + 1);

		MPI_Bcast(&Fx[0][0][0], len, MPI_MYFLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&Fy[0][0][0], len, MPI_MYFLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&Fz[0][0][0], len, MPI_MYFLOAT, MASTER, MPI_COMM_WORLD);
		MPI_Bcast(&Fp[0][0][0], len, MPI_MYFLOAT, MASTER, MPI_COMM_WORLD);
		
		goto skip_computation;
	}

	compute_ewald_correction_table();

	write_ewald_correction_table();

	skip_computation:;

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

	rprintf("done \n\n");

	return ;
}

static void compute_ewald_correction_table()
{
	rprintf("\n   Computing tables ");

	int size = p3(N_EWALD + 1) / Sim.NRank; // generic size
	int beg = Task.Rank * size;

	int len = size; // specific size treating last one

	if (Task.Rank == (Sim.NRank - 1))
		len = p3(N_EWALD + 1) - beg;

	int end = beg + len;

	int cnt = 0;

	#pragma omp parallel for
	for (int i = 0; i < N_EWALD+1; i++) {

		for (int j = 0; j < N_EWALD+1; j++) {

			for (int k = 0; k < N_EWALD+1; k++) {

				if (((cnt++) % (size/10/Sim.NThreads)) == 0)
					rprintf(".");

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

	for (int rank = 0; rank < Sim.NRank; rank++) {  // merge

		beg = rank * size;

		len = size;

		if (Task.Rank == (Sim.NRank - 1))
			len =  - beg;

		MPI_Bcast(&Fx[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fy[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fz[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fp[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);

	} // for rank

	return ;
}


static bool read_ewald_correction_table()
{
	FILE *fp = fopen(fname, "r");

	if (fp == NULL)
		return false;

	rprintf("\n   Reading Ewald tables from %s \n", fname);

	size_t nFloat = p3(N_EWALD+1);

	Fread(&Fx[0][0][0], sizeof(Fx[0][0][0]), nFloat, fp);
	Fread(&Fy[0][0][0], sizeof(Fy[0][0][0]), nFloat, fp);
	Fread(&Fz[0][0][0], sizeof(Fz[0][0][0]), nFloat, fp);

	Fread(&Fp[0][0][0], sizeof(Fp[0][0][0]), nFloat, fp);

	fclose(fp);

	return true;
}

static void write_ewald_correction_table()
{
	FILE *fp = fopen(fname, "w");

	rprintf("\n   Writing Ewald tables to %s \n", fname);

	size_t nFloat = p3(N_EWALD+1);

	Fwrite(&Fx[0][0][0], sizeof(Fx[0][0][0]), nFloat, fp);
	Fwrite(&Fy[0][0][0], sizeof(Fy[0][0][0]), nFloat, fp);
	Fwrite(&Fz[0][0][0], sizeof(Fz[0][0][0]), nFloat, fp);

	Fwrite(&Fp[0][0][0], sizeof(Fp[0][0][0]), nFloat, fp);

	fclose(fp);

	return ;
}



/*
 * Compute force correction term from Hernquist+ 1991 (2.14b, 2.16), this
 * is the difference between the force from the infinite lattice and the
 * nearest image. "x" and "force" are in boxsize units.
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

				double val = erfc(Alpha * r) + 2*Alpha * r/sqrt(PI)
							* exp(-p2(Alpha*r));

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

				if (h2 > 0) {

					double val = 2.0/h2 * exp(-PI*PI*h2/p2(Alpha))
									* sin(2*PI*hdotx);

					force[0] -= h[0] * val;
					force[1] -= h[1] * val;
					force[2] -= h[2] * val;
				}
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

				sum1 += erfc(Alpha * dr)/dr;
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

				if (h2 > 0) 
					sum2 += 1.0/(PI*h2) * exp(-PI*PI*h2 / p2(Alpha)) 
						* cos(2*PI*hdotr);
			}
		}
	}

	return PI/p2(Alpha) - sum1 - sum2 + 1.0/ALENGTH3(r);
}

#endif // GRAVITY && PERIODIC
