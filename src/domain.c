#include "globals.h"
#include "domain.h"
#include "peano.h"

#define NBUNCHES_PER_THREAD 64

static void find_global_domain();
static void fill_bunchlist();
static double compute_imbalance();
static void split_bunches();
static void update_bunchlist();
static void communicate_particles();

struct Bunch_Info *B = NULL;

/* 
 * This function distributes particles in bunches, which are continuous
 * on the Peano curve. The bunches correspond to nodes of the tree. 
 * Bunches also give the ghost nodes in the tree and are later
 * used to construct the top node tree. 
 * We keep a global list of all the bunches that also contains the 
 * workload and memory footprint of each bunch. An optimal way of 
 * distributing bunches minimises memory and workload imbalance over all
 * Tasks. 
 * */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

#ifndef PERIODIC
	find_global_domain();
#endif
 	
	Sort_Particles_By_Peano_Key();

	/*
	 while ((max_mem_imbal > 0.05) && (max_cpu_imbal > 0.05)) {
	
	  distribute();
	  decompose();

	  double max_mem_imbal = 0; 
	  double max_cpu_imbal = 0; 
	  measure(&max_mem_imbal, &max_cpu_imbal);

		
	  } 
	
	communicate_particles();
	*/

 	Sort_Particles_By_Peano_Key();

	Profile("Domain Decomposition");

	return ;
}

void Init_Domain_Decomposition()
{
/*	NBunches = fmin(floor(log2(Sim.NTask)/log2(8)), 1) * NBUNCHES_PER_THREAD;

	double nSide = log2(NBunches)/log2(8);
	
	size_t nBytes = NBunches * sizeof(B); 

	B = Malloc(nBytes, "Bunchlist");
	Tree2Bunch = Malloc(nBytes, "Tree2Bunch");
	Bunch2Tree = Malloc(nBytes, "Bunch2Tree");

	memset(B, 0, nBytes);

	const double size = 1 / nSide;

	int b = 0;

	#pragma omp parallel
	{
	
	#pragma omp single nowait
	for (int i = 0; i < nSide; i++) {

		for (int j = 0; j < nSide; j++) {
		
			for (int k = 0; k < nSide; k++) {
				
				double x =  i * size;
				double y =  j * size;
				double z =  k * size;
				
				Bunchlist[b].Key = Peano_Key(x,y,z) >> 64;

				b++;
			}
		}
	}

	find_global_domain();
	
	} // omp parallel

	rprintf("Domain Decomposition: %d Bunches @ Level %d\n"
			"   Initial size %g, origin at %4g %4g %4g \n\n", NBunches, nSide,
			Domain.Size, Domain.Origin[0],Domain.Origin[1], Domain.Origin[2]);
*/
	return;
}

static void fill_bunchlist()
{
	/*int i = 0;

	memset(&B[i].Target, 0, sizeof(*B)-sizeof(uint64_t));

	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		Float px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		Float py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		Float pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		uint64_t key = Peano_Key(px, py, pz) >> 64;

		while (B[i].Key < key) // particles are ordered
			memset(&B[i++].Target, 0, sizeof(*B)-sizeof(64/CHAR_BIT));

		B[i].Npart++;
		//B[i].CPU_Cost += P[ipart].Cost;
	
	}
*/
	return ;
}

static double compute_imbalance()
{
	return 1;
}

static void communicate_particles()
{
	return ;
}

static void update_bunchlist()
{
	return ;
}



/*
 * This finds the global domain origin and the maximum extend
 */

double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX, 
	   min_x = DBL_MIN, min_y = DBL_MIN, min_z = DBL_MIN;

static void find_global_domain()
{

	#pragma omp for reduction(max:max_x,max_y,max_z) \
		reduction(min:min_x,min_y,min_z)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		max_x = fmax(max_x, P[ipart].Pos[0]);
		max_y = fmax(max_y, P[ipart].Pos[1]);
		max_z = fmax(max_z, P[ipart].Pos[2]);

		min_x = fmin(min_x, P[ipart].Pos[0]);
		min_y = fmin(min_y, P[ipart].Pos[1]);
		min_z = fmin(min_z, P[ipart].Pos[2]);
	}
	
	#pragma omp single // do an MPI reduction
	{

	double local_max[3] = { max_x, max_y, max_z  };
	double local_min[3] = { min_x, min_y, min_z  };

	double global_max[3] = { 0 };
	double global_min[3] = { 0 };

	MPI_Allreduce(&local_max, &global_max, 3, MPI_DOUBLE, MPI_MAX,
			MPI_COMM_WORLD);

	MPI_Allreduce(&local_min, &global_min, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	Domain.Size = fabs(global_max[0] - global_min[0]);
	Domain.Size = fmax(Domain.Size, fabs(global_max[1] - global_min[1]));
	Domain.Size = fmax(Domain.Size, fabs(global_max[2] - global_min[2]));

#ifndef PERIODIC
	Sim.Boxsize[0] = Domain.Size;
	Sim.Boxsize[1] = Domain.Size;
	Sim.Boxsize[2] = Domain.Size;
#endif // PERIODIC

	Domain.Origin[0] = global_min[0];
	Domain.Origin[1] = global_min[1];
	Domain.Origin[2] = global_min[2];

	Domain.Center[0] = Domain.Origin[0] + 0.5 * Domain.Size;
	Domain.Center[1] = Domain.Origin[1] + 0.5 * Domain.Size;
	Domain.Center[2] = Domain.Origin[2] + 0.5 * Domain.Size;
	
	} // omp single

#ifdef DEBUG
	printf("\nDomain size %g, origin at %4g %4g %4g \n\n", 
			Domain.Size,
			Domain.Origin[0],Domain.Origin[1],Domain.Origin[2]);
#endif

	return ;
}
