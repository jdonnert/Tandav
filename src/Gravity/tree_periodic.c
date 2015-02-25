#include "../globals.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

#include "gravity.h"
#include "../domain.h"

static inline Float nearest(const Float dx);
static void compute_ewald_force(const int i, const int j, const int k,
								const Float r[3], double force[3]);
static double compute_ewald_potential(const Float r[3]);


const static int N_EWALD = 64; // Hernquist+ 1991
const static double alpha = 2.0;


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

	/* this follows closely Gadget-2 */

	int size = p3(N_EWALD + 1) / Sim.NRank; // generic size
	int beg = Task.Rank * size

	int len = size; // specific size treating last one

	if (Task.Rank == (Sim.NRank - 1))
		len = p3(N_EWALD + 1) - beg;

	int end = beg + len

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

				potcorr[i][j][k] = compute_ewald_potential(r);
			} // k
		} // j
	} // i

	#pragma omp parallel for schedule(static:1)
	for (int rank = 0; rank < Sim.NRank; rank+) {  // merge

		beg = rank * size;

		len = size;

		if (Task.Rank == (Sim.NRank - 1))
			len = p3(N_EWALD + 1) - beg;

		MPI_Bcast(&Fx[0][0][beg], len, MPI_MYFLOAT, task, MPI_COMM_WORLD);
		MPI_Bcast(&Fy[0][0][beg], len, MPI_MYFLOAT, task, MPI_COMM_WORLD);
		MPI_Bcast(&Fz[0][0][beg], len, MPI_MYFLOAT, task, MPI_COMM_WORLD);
		MPI_Bcast(&Fp[0][0][beg], len, MPI_MYFLOAT, task, MPI_COMM_WORLD);
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

/*
 * Compute force correction term from Hernquist+ 1992 (2.14b, 2.16), this
 * is the difference between the force from the infinite lattice and the
 * nearest image. "r" and "force" are in boxsize units.
 */

static void compute_ewald_force(const int i, const int j, const int k,
								const Float r[3], double force[3])
{
	memset(&force[0], 0, 3*sizeof(*force));

	if (i == 0 && j == 0 && k == 0)
		return ;

	double r2 = ASCALPROD3(r);

	force[0] = r[0] / (r2 * sqrt(r2));
	force[1] = r[1] / (r2 * sqrt(r2));
	force[2] = r[2] / (r2 * sqrt(r2));

	int n[3] = { 0 };

	for (n[0] = -4; n[0] < 5; n[0]++) {
		for (n[1] = -4; n[1] < 5; n[1]++) {
			for (n[2] = -4; n[2] < 5; n[2]++) {

				double dx[3] = { r[0]-n[0], r[1]-n[1], r[2]-n[2] };

				double r = ALENGTH3(dx);

				double val = erfc(alpha * r) + 2*alpha * r/sqrt(pi)
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

				double hdotr = r[0]*h[0] + r[1]*h[1] + r[2]*h[2];

				int h2 = ASCALPROD3(h);

				if (h2 <= 0)
					continue;

				double val = 2.0/h2 * exp(-pi*pi*h2/p2(alpha))
									* sin(2*pi*hdotr);

				force[0] -= h[0] * val;
				force[1] -= h[1] * val;
				force[2] -= h[2] * val;
			}
		}
	}

	return ;
}

static double compute_ewald_potential(const Float r[3])
{
	if ((r[0] == 0) && (r[1] == 0) && (r[2] == 0))
		return 2.8372975; // == U, eq. 2.15

	double sum 1 = 0;
	int n[3] = { 0 };

	for (n[0] = -4; n[0] < 5; n[0]++) {
		for (n[1] = -4; n[1] < 5; n[1]++) {
			for (n[2] = -4; n[2] < 5; n[2]++) {

				double dx[3] = { r[0]-n[0], r[1]-n[1], r[2]-n[2] };

				double r = ALENGTH3(dx);

				sum1 += erfc(alpha * r)/r;
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

				double val = 2.0 / (h2 * exp(-pi*pi*h2/p2(alpha))
									* sin(2*pi*hdotr));

				sum2 += 1.0/(pi*h2) * exp(-pi*pi*h2 / p2(alpha)) 
					* cos(2*pi*hdotr);

			}
		}
	}

	double r = ALENGTH3(r);

	return pi / p2(alpha) - sum1 - sum2 + 1.0/r;
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



#endif
