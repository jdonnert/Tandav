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

static void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--) {
		
		printf("%llu", (val & (1ULL << i) ) >> i);
		
		if (i % 3 == 0 && i != 0)
			printf(".");
	}
	printf("\n");fflush(stdout);

	return ;
}

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

	find_global_domain();
 	
	Sort_Particles_By_Peano_Key();

	fill_bunchlist(0, 0, NBunches, Task.Npart_Total);

	for (int i = 0; i < NBunches; i++) {
		printf("%d %d \n", i, B[i].Npart);
		print_int_bits64(B[i].Key);
	}
	exit(0);
	

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
	find_global_domain();

	NBunches = 8;

	int nSide = log2(NBunches)/log2(8);
	
	size_t nBytes = NBunches * sizeof(*B); 

	B = Malloc(nBytes, "Bunchlist");

	memset(B, 0, nBytes);

	const int shift = 63 - log2(NBunches);

	#pragma omp parallel for 
	for (int i = 0; i < NBunches; i++) 
		B[i].Key = (uint64_t) (i+1) << shift;  // all bit permutations 

	rprintf("Domain Decomposition: \n"
			"   %d Bunches @ Level %d for %d Tasks\n"
			"   Initial size %g, origin at %4g %4g %4g \n\n", 
			NBunches, nSide, Sim.NTask, Domain.Size, Domain.Origin[0],
			Domain.Origin[1], Domain.Origin[2]);

	return;
}

static void fill_bunchlist(const int first_bunch, const int first_part, 
		const int nBunch, const nPart)
{
	const int last_bunch = first_bunch + nBunch + 1;
	const int last_part = first_part + nPart + 1;

	#pragma omp for
	for (int i = first_bunch; i < last_bunch; i++)
		B[i].Target =  B[i].Npart = B[i].CPU_Cost = 0;

printf("DING ! \n");	
	int i = first_bunch;

	#pragma omp for
	for (int ipart = first_part; ipart < last_part; ipart++) {
		
		Float px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		Float py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		Float pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		uint64_t pkey = Peano_Key(px, py, pz) >> 64;

		while ((B[i].Key < pkey) && (i < last_bunch))  // particles are ordered
			i++;

		printf("%d %d %d ", ipart, i , B[i].Npart);
		print_int_bits64(B[i].Key); printf("\n");

		#pragma omp critical 
		{

		if (B[i].Npart == 0) {
		
			B[i].Pos[0] = P[ipart].Pos[0];
			B[i].Pos[1] = P[ipart].Pos[1];
			B[i].Pos[2] = P[ipart].Pos[2];
		}

		B[i].Npart++;

		//B[i].CPU_Cost += P[ipart].Cost;
		} // omp critical

	}

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
 * Find the global domain origin and the maximum extent
 */

double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX, 
	   min_x = DBL_MIN, min_y = DBL_MIN, min_z = DBL_MIN;

static void find_global_domain()
{

#ifdef PERIODIC
	
	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;

	Domain.Size = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));
	
	Domain.Center[0] = Domain.Origin[0] + Domain.Size/2;
	Domain.Center[1] = Domain.Origin[1] + Domain.Size/2;
	Domain.Center[2] = Domain.Origin[2] + Domain.Size/2;

#else // ! PERIODIC

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

	Sim.Boxsize[0] = Sim.Boxsize[1] = Sim.Boxsize[2] = Domain.Size;

	Domain.Origin[0] = global_min[0];
	Domain.Origin[1] = global_min[1];
	Domain.Origin[2] = global_min[2];

	Domain.Center[0] = Domain.Origin[0] + 0.5 * Domain.Size;
	Domain.Center[1] = Domain.Origin[1] + 0.5 * Domain.Size;
	Domain.Center[2] = Domain.Origin[2] + 0.5 * Domain.Size;
	
	} // omp single

#endif // PERIODIC

#ifdef DEBUG
	printf("\nDomain size %g, origin at %4g %4g %4g \n\n", 
			Domain.Size,
			Domain.Origin[0],Domain.Origin[1],Domain.Origin[2]);
#endif

	return ;
}
