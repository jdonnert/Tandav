#include "globals.h"
#include "domain.h"
#include "peano.h"

static void find_global_domain();
static void fill_bunchlist();
static double compute_imbalance();
static void split_bunches();
static void update_bunchlist();
static void communicate_particles();

/* This function distributes particles in bunches, which are continuous
 * on the Peano curve. This way, the bunches correspond to parts
 * of the tree. Bunches also give the ghost nodes in the tree and are later
 * used to construct the top node tree. 
 * We keep a global list of all the bunches that also contains the workload and 
 * memory footprint of each bunch. An optimal way of distributing bunches
 * minimises memory and workload imbalance over all Tasks. */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

	find_global_domain();

 	Sort_Particles_By_Peano_Key();

	/*for (;;) {
	
		fill_bunchlist();

		double imbalance = compute_imbalance();

		if ( imbalance < 0.05 )
			break;
		
		split_bunches();
	} 
	
	communicate_particles();
*/
	Profile("Domain Decomposition");

	return ;
}

void Init_Domain_Decomposition()
{
	#pragma omp parallel
	find_global_domain();

	/*int nTask = Sim.NTask, nThreads = Sim.NThreads;
	
	size_t nBytes = (1ULL << DOMAIN_CEILING) * Sim.NRank * sizeof(Bunchlist); 

	Bunchlist = Malloc(nBytes);

	memset(Bunchlist, 0, nBytes);

	const uint32_t nBunches = 1UL << (int) ceil(log2(nTask));
	const uint32_t nSide = log2(nBunches);
	const double size = Sim.Boxsize[0] / nSide;

	for (int i = 0; i < nSide; i++){

		for (int j = 0; j < nSide; j++){
		
			for (int k = 0; k < nSide; k++) {
				
				size_t idx = i * p2(nSide) + j * nSide + k;
				
				Bunchlist[idx].Pos[0] = i * size;
				Bunchlist[idx].Pos[1] = i * size;
				Bunchlist[idx].Pos[2] = i * size;
			}
		}
	}

	for (int i = 0; i < nBunches; i++) {
	
		Float x = Bunchlist[i].Pos[0];
		Float y = Bunchlist[i].Pos[1];
		Float z = Bunchlist[i].Pos[2];
		
		Bunchlist[i].Key = Peano_Key(x,y,z, Sim.Boxsize);
	}
*/
	return;
}
static void fill_bunchlist()
{

	return ;
}

static double compute_imbalance()
{
	return 1;
}

static void split_bunches()
{

	return ;
}

static void communicate_particles()
{
	return ;
}

static void update_bunchlist()
{
	return ;
}

static double global_min[3] = { 0 }; 
static double global_max[3] = { 0 };

static void find_global_domain()
{
	double local_max[3] = { 0 };
	double local_min[3] = { 0 };

	#pragma omp for nowait
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		local_max[0] = fmax(local_max[0], P[ipart].Pos[0]);
		local_max[1] = fmax(local_max[1], P[ipart].Pos[1]);
		local_max[2] = fmax(local_max[2], P[ipart].Pos[2]);

		local_min[0] = fmin(local_min[0], P[ipart].Pos[0]);
		local_min[1] = fmin(local_min[1], P[ipart].Pos[1]);
		local_min[2] = fmin(local_min[2], P[ipart].Pos[2]);
	}

	#pragma omp critical // do an omp reduction
	{

	global_max[0] = fmax(global_max[0], local_max[0]);
	global_max[1] = fmax(global_max[1], local_max[1]);
	global_max[2] = fmax(global_max[2], local_max[2]);
	
	global_min[0] = fmin(global_min[0], local_min[0]);
	global_min[1] = fmin(global_min[1], local_min[1]);
	global_min[2] = fmin(global_min[2], local_min[2]);
	
	} // omp critical

	#pragma omp single // do an MPI reduction
	{

	memcpy(local_max, global_max, 3*sizeof(*local_max)); // copy back as buf
	memcpy(local_min, global_min, 3*sizeof(*local_min));

	MPI_Allreduce(&local_max, &global_max, 3, MPI_DOUBLE, MPI_MAX,
			MPI_COMM_WORLD);

	MPI_Allreduce(&local_min, &global_min, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	Domain.Size[0] = fabs(global_max[0] - global_min[0]);
	Domain.Size[1] = fabs(global_max[1] - global_min[1]);
	Domain.Size[2] = fabs(global_max[2] - global_min[2]);

#ifndef PERIODIC
	Sim.Boxsize[0] = Domain.Size[0];
	Sim.Boxsize[1] = Domain.Size[1];
	Sim.Boxsize[2] = Domain.Size[2];
#endif // PERIODIC

	Domain.Corner[0] = global_min[0];
	Domain.Corner[1] = global_min[1];
	Domain.Corner[2] = global_min[2];

#ifdef DEBUG
	printf("\n%d, %d found domain size %g %g %g at %g %g %g \n\n", 
			Task.Rank, Task.Thread_ID,
			Domain.Size[0], Domain.Size[1],Domain.Size[2],
			Domain.Corner[0],Domain.Corner[1],Domain.Corner[2]);
#endif

	} // omp single

	return ;
}
