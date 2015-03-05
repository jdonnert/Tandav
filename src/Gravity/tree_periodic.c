#include "../globals.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

#include "gravity.h"
#include "../domain.h"

#define N_EWALD 64 // Hernquist+ 1991

static void ewald_correction(const Float dr[3], Float f[3]);
static void compute_ewald_correction_table();
static void write_ewald_correction_table();
static bool read_ewald_correction_table();
static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
										Float* Accel, Float *Pot);
static void compute_ewald_force(const int, const int, const int,
								const double r[3], double force[3]);
#ifdef OUTPUT_GRAV_POTENTIAL    
static double compute_ewald_potential(const double r[3]);
static void ewald_potential(const Float dr[3], Float p[1]);
#else
static inline void compute_ewald_potential(const double r[3]) {};
static inline void ewald_potential(const Float dr[3], Float p[1]) {};
#endif // OUTPUT_GRAV_POTENTIAL

const static double Alpha = 2.0; // Hernquist+ 1992 (2.13)

static double Box2Ewald_Grid = 0;
static double Boxsize = 0;
static double Boxhalf = 0;

static Float Fx[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fy[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fz[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };
static Float Fp[N_EWALD+1][N_EWALD+1][N_EWALD+1] = { 0 };


/*
 * Compute the correction to the gravitational force due the periodic
 * infinite box using the tree and the Ewald method (Hernquist+ 1992).
 */

void Gravity_Tree_Periodic()
{
	Profile("Grav Tree Periodic");

	Profile("Grav Tree Periodic");
	return ;
}

/*
 * This walks a subtree starting at "tree_start" to estimate the Ewald 
 * correction to the gravitational force. Because the force is proportional to
 * the distance a different opening criterion can be used, similar to Gadget-2.
 * However, this cannot happen across a periodic boundary. 
 */

static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
		Float Accel[3], Float Pot[1])
{
	const Float fac = ALENGTH3(P[ipart].Acc) / Const.Gravity
		* TREE_OPEN_PARAM_REL;
	
	const Float pos_i[3] = {P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2]};
	
	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		bool want_open_node = false;

		if (Tree[node].DNext < 0) { // particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++) {

				if (jpart == ipart)
					continue;

				Float dr[3] = { pos_i[0] - P[jpart].Pos[0],
							    pos_i[1] - P[jpart].Pos[1],
					            pos_i[2] - P[jpart].Pos[2] };

				Tree_Periodic_Nearest(dr);

				ewald_correction(dr, Accel); 

 				ewald_potential(dr, Pot);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Tree_Periodic_Nearest(dr);

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = node_size(node); // now check opening criteria

		if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion
			
			want_open_node = true;

			goto skip;
		}
		
		dr[3] = { pos_i[0] - Tree[node].Pos[0],
			      pos_i[1] - Tree[node].Pos[1],
			      pos_i[2] - Tree[node].Pos[2] };

		if (dr[0] < 0.6 * nSize) { // particle inside node ? 
			
			if (dr[1] < 0.6 * nSize) {

				if (dr[2] < 0.6 * nSize) {
				
					want_open_node = true;
				
				}
			}
		}

		skip:;

		if (want_open_node) { // check if we really have to open

			Tree_Periodic_Nearest(dr);

			Float size = Node_Size(node);
			
			if (fabs(dr[0]) > 0.5 * (Boxsize - size)) { 
				
				node++;

				continue;
			}

			if (fabs(dy) > 0.5 * (Boxsize - size)) {
				
				node++;

				continue;
			
			}

			if (fabs(dz) > 0.5 * (Boxsize - size)) {
				
				node++;

				continue;
			
			}

			if (Tree[node].Size > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}
		
		} // if want_open_node

		ewald_correction(dr, accel); // use node

 		ewald_potential(dr, pot); // OUTPUT_GRAVITATIONAL_POTENTIAL

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}



void Gravity_Tree_Periodic_Init()
{
	rprintf("Init Ewald correction ");

	Boxsize = Sim.Boxsize[0];
	Boxhalf = Boxsize / 2;

	Box2Ewald_Grid = 2 * N_EWALD / Boxsize;

	bool table_found = false;

	if (Task.Is_Master)
		 table_found = read_ewald_correction_table();

	MPI_Bcast(&table_found, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	if (table_found)
		goto skip_computation;

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

/*
 * Do the periodic mapping on a 3D distance array. We rely on link time 
 * optimization of the compiler to do the inlining for us. Make sure to put
 * the appropriate compiler switches.
 */

void Tree_Periodic_Nearest(Float dr[3])
{
	for (int i = 0; i < 3; i++) {
	
		if (dr[i] > Boxhalf)
			dr[i] -= Boxsize;
		else if (dr[i] < -Boxhalf)
			dr[i] += Boxsize;
	}

	return ;
}

/* 
 * this follows closely Gadget-2 (Springel 2006)
 */

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
			len = p3(N_EWALD + 1) - beg;

		MPI_Bcast(&Fx[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fy[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fz[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);
		MPI_Bcast(&Fp[0][0][beg], len, MPI_MYFLOAT, rank, MPI_COMM_WORLD);

	} // for rank

	return ;
}

/*
 * read/write the ewald tables from/to a binary file
 */

const static char fname[CHARBUFSIZE] = { "ewald_table.dat" };

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

	rprintf("\n   Writing Ewald tables from %s \n", fname);

	size_t nFloat = p3(N_EWALD+1);

	Fwrite(&Fx[0][0][0], sizeof(Fx[0][0][0]), nFloat, fp);
	Fwrite(&Fy[0][0][0], sizeof(Fy[0][0][0]), nFloat, fp);
	Fwrite(&Fz[0][0][0], sizeof(Fz[0][0][0]), nFloat, fp);

	Fwrite(&Fp[0][0][0], sizeof(Fp[0][0][0]), nFloat, fp);

	fclose(fp);

	return ;
}

/*
 * Get Ewald correction from grid using CIC binning (Hockney & Eastwood)
 */

static void ewald_correction(const double dr[3], double f[3])
{
	f[0] = f[1] = f[2] = 0;

	double dx = dr[0];
	int signx = -1;

	if (dx < 0) {
	
		dx = -dx;
		signx = 1;
	}

	double dy = dr[1];
	int signy = -1;

	if (dy < 0) {
	
		dy = -dy;
		signy = 1;
	}

	double dz = dr[2];
	int signz = -1;

	if (dz < 0) {
	
		dz = -dz;
		signz = 1;
	}

	double u = dx * Box2Ewald_Grid;
	double v = dy * Box2Ewald_Grid;
	double w = dz * Box2Ewald_Grid;

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

#ifdef GRAVITY_POTENTIAL
static void ewald_potential(const double dr[3], double p[1])
{
	p[0] = 0;

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
#endif // PERIODIC

/*
 * Compute force correction term from Hernquist+ 1992 (2.14b, 2.16), this
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

				if (h2 <= 0)
					continue;

				double val = 2.0/h2 * exp(-PI*PI*h2/p2(Alpha))
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

				if (h2 <= 0)
					continue;

				sum2 += 1.0/(PI*h2) * exp(-PI*PI*h2 / p2(Alpha)) 
					* cos(2*PI*hdotr);
			}
		}
	}

	return PI/p2(Alpha) - sum1 - sum2 + 1.0/ALENGTH3(r);
}

#undef N_EWALD

#endif
