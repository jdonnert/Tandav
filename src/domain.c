#include "globals.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "peano.h"

#define DOMAIN_SPLIT_MEM_THRES (0.0)
#define DOMAIN_NBUNCHES_PER_THREAD 16

static void find_global_domain();
static void fill_bunches(const int, const int, const int, const int);
static void remove_empty_bunches();
static bool bunch_is_overloaded(const int b);
static int split_bunch(const int i);
static void find_imbalance();
static void communicate_particles();
static void communicate_bunch_list();

static int compare_bunches_by_key(const void *a, const void *b);
static int compare_bunches_by_target(const void *a, const void *b); 

union Domain_Node_List *D; 

static double max_mem_imbal = DBL_MAX;
static double max_cpu_imbal = DBL_MAX;

/* 
 * Distribute particles in bunches, which are continuous
 * on the Peano curve. The bunches correspond to nodes of the tree. 
 * Bunches also give the ghost nodes in the tree. 
 * We keep a global list of all the bunches that also contains the 
 * workload and memory footprint of each bunch. An optimal way of 
 * distributing bunches minimises memory and workload imbalance over all
 * Tasks. 
 */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

	find_global_domain();
 	
	Sort_Particles_By_Peano_Key();

	#pragma omp for
	for (int i = 0; i < NBunches; i++) { // reset bunchlist

		D[i].Bunch.Target =  D[i].Bunch.Npart = D[i].Bunch.Cost = 0;
		D[i].Bunch.First_Part = INT_MAX;
	}

	// restore_complete_bunchlist();

	fill_bunches(0, NBunches, 0, Task.Npart_Total); // let's see what we have

	for (int i = 0; i < NBunches; i++) 
		printf("%d %d %d \n", i, D[i].Bunch.Npart, D[i].Bunch.First_Part);

printf("-----------------------\n");

	remove_empty_bunches();

	Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

for (int i = 0; i < NBunches; i++) {
printf("%d | %5d | ", i, D[i].Bunch.Npart);
Print_Int_Bits64(D[i].Bunch.Key);
printf("            ");
int ifirst = D[i].Bunch.First_Part;

Float px = (P[ifirst].Pos[0] - Domain.Origin[0]) / Domain.Size;
Float py = (P[ifirst].Pos[1] - Domain.Origin[1]) / Domain.Size;
Float pz = (P[ifirst].Pos[2] - Domain.Origin[2]) / Domain.Size;
peanoKey pkey = Peano_Key(px, py, pz);
Print_Int_Bits128(pkey);
		
int ilast = ifirst + D[i].Bunch.Npart - 1;
printf("            ");
 px = (P[ilast].Pos[0] - Domain.Origin[0]) / Domain.Size;
 py = (P[ilast].Pos[1] - Domain.Origin[1]) / Domain.Size;
 pz = (P[ilast].Pos[2] - Domain.Origin[2]) / Domain.Size;
 pkey = Peano_Key(px, py, pz);
Print_Int_Bits128(pkey);
}
printf("-----------------------\n");

	int max_level = 1;

	while ((max_mem_imbal > 0.05) || (max_cpu_imbal > 0.05)) {
	
		int old_nBunches = NBunches;
	
		for (int i = 0; i < old_nBunches; i++ ) {
	
			if (bunch_is_overloaded(i)) { // split into 8
	printf("Split %d \n", i);
				int first_new_bunch = split_bunch(i);

				fill_bunches(first_new_bunch, 8, D[i].Bunch.First_Part, 
						D[i].Bunch.Npart);

				max_level = fmax(max_level, 1 + D[i].Bunch.Level);

				memset(&D[i].Bunch, 0, sizeof(*D));
			}
		}	
		
		// communicate_bunches();

		remove_empty_bunches();

	//Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_target);
		// distribute();

		Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

		find_imbalance();

		if (max_level == N_SHORT_TRIPLETS-1)
			break;

	} // while
	
	
int sum = 0;
for (int i = 0; i < NBunches; i++) {
	sum+= D[i].Bunch.Npart;
printf("%d | %5d %5d %d  ", i, D[i].Bunch.Npart, sum, D[i].Bunch.Level);
Print_Int_Bits64(D[i].Bunch.Key);
}
printf("\n-----------------------\n");

	Profile("Domain Decomposition");

	return ;
}

void Init_Domain_Decomposition()
{
	find_global_domain();

	NBunches = 1;

	size_t nBytes = 8 * Sim.NTask* DOMAIN_NBUNCHES_PER_THREAD * sizeof(*D); 

	D = Malloc(nBytes, "Bunch Tree List");

	memset(D, 0, nBytes);

	D[0].Bunch.Key = 0xFFFFFFFFFFFFFFFF;
	D[0].Bunch.Npart = Sim.Npart_Total;
	D[0].Bunch.Level = D[0].Bunch.Target = 0;

	int new = split_bunch(0); // make first 8, new = 1

	memmove(&D[0], &D[new], 8*sizeof(*D));
	NBunches--;

	Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

	rprintf("Domain Decomposition: \n"
			"   Initial size %g, origin at %4g %4g %4g \n\n", 
			Domain.Size, Domain.Origin[0], Domain.Origin[1], Domain.Origin[2]);

	return;
}

/*
 * Transform the top nodes back into a bunch list. Add nodes so the complete 
 * domain is covered.
 */

void restore_complete_bunchlist()
{

	return ;
}

/*
 * We split a bunch into 8 sub-bunches/nodes, adding the largest peano key 
 * contained in the bunch. The position of the bunch is set during filling to 
 * a random particle position. From it we can later construct the top node 
 * position during Tree construction.
 */

static int split_bunch(const int parent)
{
	int first = 0; 

	#pragma omp critical // add new nodes at the end
	{

	first = NBunches; 

	NBunches += 8;

	} // omp critical

	for (int i = 0; i < 8; i++) {

		int dest = first + i;

		D[dest].Bunch.Level = D[parent].Bunch.Level + 1;

		int shift = N_SHORT_BITS - 1 - 3 * D[dest].Bunch.Level;

		shortKey bitmask = (1ULL << (N_SHORT_BITS - 1)) | (0x7ULL << shift);
		shortKey keyfragment = ((shortKey) i) << shift;

		D[dest].Bunch.Key = (D[parent].Bunch.Key & ~bitmask) | keyfragment;
		D[dest].Bunch.Npart = 0;
		D[dest].Bunch.First_Part = INT_MAX;
	}

	return first;
}

static void remove_empty_bunches()
{
	#pragma omp single
	{

	int old_nBunches = NBunches;

	int i = 0;

	while (i < NBunches) {
		
		if (D[i].Bunch.Npart == 0) {  // remove
	
			memmove(&D[i], &D[i+1], (NBunches - i - 1)*sizeof(*D));

			NBunches--;

			continue;
		} 

		i++;
	
	}

	printf("(%d:%d) Removed %d bunches, now %d Bunches\n", Task.Rank, 
			Task.Thread_ID, old_nBunches-NBunches, NBunches);

	} // omp single

	return ;
}

/*
 * distribute local particles into bunches
 */

static void fill_bunches(const int first_bunch, const int nBunch, 
		 const int first_part, const int nPart)
{
	const int last_part = first_part + nPart;
	
	int i = first_bunch;

	int npart = 0;
	float cost = 0;
	int first = INT_MAX;

	#pragma omp for
	for (int ipart = first_part; ipart < last_part; ipart++) {
		
		double px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		double py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		double pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		shortKey pkey = Short_Peano_Key(px, py, pz);

		while (D[i].Bunch.Key < pkey) { // particles are ordered by key
			
			if (npart > 0) {
				
				#pragma omp atomic update
				D[i].Bunch.Npart += npart;

				#pragma omp atomic update
				D[i].Bunch.Cost += cost;

				#pragma omp critical
				if (first < D[i].Bunch.First_Part)
					D[i].Bunch.First_Part = first;
			
				npart = cost = 0;
				first = INT_MAX;
			}

			i++;
		}

		npart++;
		cost += P[ipart].Cost;
		first = fmin(ipart, first);
	}

	if (npart > 0) { // dont forget the last bunch !
		
		#pragma omp atomic update
		D[i].Bunch.Npart += npart;

		#pragma omp atomic update
		D[i].Bunch.Cost += cost;

		#pragma omp critical
		if (first < D[i].Bunch.First_Part)
			D[i].Bunch.First_Part = first;
	}
	
	return ;
}



/*
 * This function defines the metric that decides if a bunch has to be refined
 * into eight sub-bunches.
 */

static void find_imbalance()
{
	max_mem_imbal = max_cpu_imbal = 0;

	const double mean_npart = Sim.Npart_Total
		/(Sim.NTask + Sim.NThreads * DOMAIN_NBUNCHES_PER_THREAD);

	#pragma omp for reduction(max:max_mem_imbal, max_cpu_imbal)
	for (int i = 0; i < NBunches; i++ ) {

		double rel_mem_load = 
			(double) (D[i].Bunch.Npart - mean_npart) / mean_npart;

		max_mem_imbal = fmax(max_mem_imbal, rel_mem_load);

		max_cpu_imbal = 0;
	}

	return ;
}

static bool bunch_is_overloaded(const int b)
{
	if (D[b].Bunch.Npart == 0)
		return false;

	const double mean_npart = Sim.Npart_Total
		/(Sim.NTask + Sim.NThreads * DOMAIN_NBUNCHES_PER_THREAD);
	
	double rel_mem_load = (double)(D[b].Bunch.Npart - mean_npart) / mean_npart;

printf("Test Bunch %d, np=%5d mean=%g delta=%g \n",
b, D[b].Bunch.Npart, mean_npart, rel_mem_load );
	
	if (rel_mem_load > DOMAIN_SPLIT_MEM_THRES)
		return true;

	return false;
}

static int compare_bunches_by_key(const void *a, const void *b) 
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Key > y->Key) - (x->Key < y->Key);
}

static int compare_bunches_by_target(const void *a, const void *b) 
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Target > y->Target) - (x->Target < y->Target);
}

static void communicate_particles()
{
	return ;
}

static void communicate_bunch_list()
{
	return ;
}
	
/*
 * Find the global domain origin, center and the maximum extent
 */

double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX, 
	   min_x = DBL_MAX, min_y = DBL_MAX, min_z = DBL_MAX;

static void find_global_domain()
{

#ifdef PERIODIC
	
	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;

	Domain.Size = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));

	for (int i = 0; i < 3; i++)
		Domain.Center[i] = Domain.Origin[i] + Domain.Size/2;

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
	
	double global_max[3] = { max_x, max_y, max_z  };
	double global_min[3] = { min_x, min_y, min_z  };

	#pragma omp single // do an MPI reduction
	{

		MPI_Allreduce(MPI_IN_PLACE, &global_max, 3, MPI_DOUBLE, MPI_MAX,
			MPI_COMM_WORLD);

		MPI_Allreduce(MPI_IN_PLACE, &global_min, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);
	
	} // omp single

	#pragma omp flush
	
	Domain.Size = fabs(global_max[0] - global_min[0]);
	Domain.Size = fmax(Domain.Size, fabs(global_max[1] - global_min[1]));
	Domain.Size = fmax(Domain.Size, fabs(global_max[2] - global_min[2]));

	for (int i = 0; i < 3; i++) {
	
		Domain.Origin[i] = global_min[i]; 
		Domain.Center[i] = Domain.Origin[i] + 0.5 * Domain.Size;
	}
#endif // PERIODIC

#ifdef DEBUG
	rprintf("\nDomain size %g, origin at %4g %4g %4g \n\n", Domain.Size,
			Domain.Origin[0], Domain.Origin[1], Domain.Origin[2]);
#endif

	return ;
}
