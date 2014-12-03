#include "globals.h"
#include "domain.h"
#include "peano.h"

#define DOMAIN_SPLIT_MEM_THRES (0.05)

static void find_global_domain();
static void fill_bunches(const int, const int, const int, const int);
static bool bunch_overloaded(const int b);
static int split_bunch(const int i);
static void communicate_particles();
static int compare_bunches(const void *a, const void *b);

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

	find_global_domain();
 	
	Sort_Particles_By_Peano_Key();

	#pragma omp for
	for (int i = 0; i < NBunches; i++) // reset bunchlist
		B[i].Target =  B[i].Npart = B[i].CPU_Cost = 0;

	fill_bunches(0, NBunches, 0, Task.Npart_Total); // let's see what we have

	for (int i = 0; i < NBunches; i++) {
		printf("%d %d \n", i, B[i].Npart);
		Print_Int_Bits64(B[i].Key);
	}
printf("-----------------------\n");
	// while (total_imbalance < X) {
	
	int old_nBunches = NBunches;
	
	for (int i = 0; i < old_nBunches; i++ ) {
	
		if (bunch_overloaded(i)) { // split into 8
printf("split %d \n", i);	
			int first_new_bunch = split_bunch(i);

			fill_bunches(first_new_bunch, 8, B[i].Target, B[i].Npart);

			memset(&B[i], 0, sizeof(*B));
		}
	}

	for (int i = 0; i < NBunches; i++) {
		printf("%d %d \n", i, B[i].Npart);
		Print_Int_Bits64(B[i].Key);
	}
	
	Qsort(Sim.NThreads, B, NBunches, sizeof(*B), &compare_bunches);

	if (old_nBunches - NBunches == 0) // reuse old distribution 
		goto skip;

	/*
	 while ((max_mem_imbal > 0.05) && (max_cpu_imbal > 0.05)) {

		
	  } 
	
	communicate_particles();
	*/

	skip:;

 	Sort_Particles_By_Peano_Key();

	Profile("Domain Decomposition");
exit(0);
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

	const double size = 0.5 * Domain.Size;

	#pragma omp parallel for 
	for (int i = 0; i < NBunches; i++) { 
	
		B[i].Level = 1;

		double px =  0.25 + 0.5 * (i & 1); // position in domain units
		double py =  0.25 + 0.5 * (i & 2) >> 1;
		double pz =  0.25 + 0.5 * (i & 4) >> 2;

		B[i].Key = Peano_Key(px, py, pz);   

		B[i].Pos[0] = Domain.Corner[0] + px * Domain.Size;
		B[i].Pos[1] = Domain.Corner[1] + py * Domain.Size;
		B[i].Pos[2] = Domain.Corner[2] + pz * Domain.Size;
		
	}

	rprintf("Domain Decomposition: \n"
			"   %d Bunches @ Level %d for %d Tasks\n"
			"   Initial size %g, origin at %4g %4g %4g \n\n", 
			NBunches, nSide, Sim.NTask, Domain.Size, Domain.Origin[0],
			Domain.Origin[1], Domain.Origin[2]);

	return;
}

static void fill_bunches(const int first_bunch, const int nBunch, 
		 const int first_part, const int nPart)
{
	const int last_part = first_part + nPart + 1;
	
	int i = first_bunch;

	#pragma omp for
	for (int ipart = first_part; ipart < last_part; ipart++) {
		
		Float px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		Float py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		Float pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		shortKey pkey = Peano_Key(px, py, pz) >> 64;

		while (B[i].Key < pkey) // particles are ordered
			i++;

		#pragma omp atomic
		B[i].Npart++;

		//#pragma omp atomic
		//B[i].CPU_Cost += P[ipart].Cost;

	}

	return ;
}

/*
 * This function defines the metric that decides if a bunch has to be refined
 * into eight sub-bunches.
 */

static bool bunch_overloaded(const int b)
{
	bool result = false;

	const int64_t mean_npart = Sim.Npart_Total / NBunches;

	double rel_mem_load = (double)(B[b].Npart - mean_npart) / mean_npart;
printf("%d %g \n",b , rel_mem_load);
	if (rel_mem_load > 1 + DOMAIN_SPLIT_MEM_THRES)
		result = true;

	return result;
}

static int split_bunch(const int src)
{
	#pragma omp atomic 
	NBunches += 8;

	const int first = NBunches - 8;

	for (int i = ; i < 8; i++) {

		int dest = first + i;

		B[dest].Size = 0.5 * B[src].Size;

		B[dest].Pos[0] = B[src].Pos[0] + (-0.5 + (i & 1)) * B[dest].Size;
		B[dest].Pos[1] = B[src].Pos[1] + (-0.5 + (i & 2) >> 1) * B[dest].Size;
		B[dest].Pos[2] = B[src].Pos[2] + (-0.5 + (i & 4) >> 2) * B[dest].Size;

		double bx = (B[dest].Pos[0] - Domain.Corner[0] ) / Domain.Size;
		double by = (B[dest].Pos[1] - Domain.Corner[1] ) / Domain.Size;
		double bz = (B[dest].Pos[2] - Domain.Corner[2] ) / Domain.Size;

		B[i].Key = Peano_Key(bx, by, bz) >> 64;

	}

	return dest;
}

static int compare_bunches(const void *a, const void *b) 
{
	const struct Bunch_Info *x = (const struct Bunch_Info *) a;
	const struct Bunch_Info *y = (const struct Bunch_Info *) b;

	return (int) (x->Key > y->Key) - (x->Key < y->Key);
}


static void communicate_particles()
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
