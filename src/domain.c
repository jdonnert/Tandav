#include "globals.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "peano.h"

#define DOMAIN_SPLIT_MEM_THRES -0.6
#define DOMAIN_NBUNCHES_PER_THREAD 4.0

static void reset_bunchlist();
static void find_global_domain();
static void fill_bunches(const int, const int, const int, const int);
static void remove_empty_bunches();
static int split_bunch(const int i);
static void reallocate_topnodes();
static void distribute();
static void communicate_particles();
static void communicate_bunches();
static void print_domain_decomposition (const int);

static int compare_bunches_by_key(const void *a, const void *b);
static int compare_bunches_by_target(const void *a, const void *b); 
static int compare_bunches_by_npart(const void *a, const void *b); 

union Domain_Node_List *D; 

static bool Keep_Splitting = false;
static double max_mem_imbal = 0, max_cpu_imbal = 0;

double Top_Node_Alloc_Factor = 0.001;
int Max_NBunches = 9;

/* 
 * Distribute particles in bunches, which are continuous
 * on the Peano curve. The bunches correspond to nodes of the tree. 
 * Bunches also give the top nodes in the tree, with some info added. 
 * We keep a global list of all the bunches that also contains the 
 * workload and memory footprint of each bunch. An optimal way of 
 * distributing bunches minimises memory and workload imbalance over all
 * Tasks. 
 * To achieve this, we measure mem & cpu cost and refine bunches until the
 * heaviest have roughly equal cost. Then we distribute them across MPI ranks
 */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

	reset_bunchlist();

	find_global_domain();
 	
	Sort_Particles_By_Peano_Key();
		
	fill_bunches(0, NBunches, 0, Task.Npart_Total); // let's see what we have

	remove_empty_bunches();

	Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

	distribute();

	int max_level = 1;

	while (Keep_Splitting && (max_level < N_SHORT_TRIPLETS-1)) {
	
		int old_nBunches = NBunches;

		for (int i = 0; i < old_nBunches; i++ ) {

			if (D[i].Bunch.Is_To_Be_Split) { // split into 8

				if (NBunches + 8 >= Max_NBunches) // make more space !
					break;

				int first_new_bunch = split_bunch(i);

				max_level = fmax(max_level, 1 + D[i].Bunch.Level);

				fill_bunches(first_new_bunch, 8, D[i].Bunch.First_Part, 
						D[i].Bunch.Npart);
				
				memset(&D[i].Bunch, 0, sizeof(*D));
			}
		}

		#pragma omp barrier
		
		if (NBunches + 8 >= Max_NBunches) { // make more space !
		
			reallocate_topnodes();

			continue;
		}
		
		communicate_bunches();

		remove_empty_bunches();

		Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

		distribute();
	} // while

	rprintf("\nDomain: %d Bunches/Top Nodes, max level %d \n\n", 
			NBunches, max_level);

#ifdef DEBUG
	print_domain_decomposition(max_level);
#endif

	communicate_particles();

	Profile("Domain Decomposition");

	return ;
}


void Init_Domain_Decomposition()
{
	find_global_domain();

	Max_NBunches = 9;

	Top_Node_Alloc_Factor = (double) Max_NBunches / Task.Npart_Total;

	size_t nBytes = Max_NBunches * sizeof(*D); 

	D = Malloc(nBytes, "D");

	memset(D, 0, nBytes);

	NBunches = 1;

	D[0].Bunch.Key = 0xFFFFFFFFFFFFFFFF;
	D[0].Bunch.Npart = Sim.Npart_Total;
	D[0].Bunch.Level = D[0].Bunch.Target = 0;
	
	int new = split_bunch(0); // make first 8, new = 1

	memmove(&D[0], &D[new], 8*sizeof(*D)); // move left by one

	NBunches--;

	Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

	rprintf("\nDomain size %g, \n   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Com    at x = %4g, y = %4g, z = %4g. \n", Domain.Size,
			Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center_Of_Mass[0], Domain.Center_Of_Mass[1],
			Domain.Center_Of_Mass[2]);

	return;
}

/*
 * This increases the room for Bunches/Topnodes by 20 %, so we can stay minimal
 */

static void reallocate_topnodes()
{
	#pragma omp single
	{
	
	Top_Node_Alloc_Factor *= 1.2;
	
	Max_NBunches = Sim.Npart_Total* Top_Node_Alloc_Factor;

	size_t nBytes = Max_NBunches * sizeof(*D); 

	mprintf("Increasing Top Node Memory to %g MB, Max %d Nodes, Factor %g \n"
			, nBytes/1024.0/1024.0, Max_NBunches, Top_Node_Alloc_Factor);

	D = Realloc(D, nBytes, "Bunchlist");

	} // omp single

	return ;
}

/*
 * Transform the top nodes back into a bunch list. Add nodes so the complete 
 * domain is covered. We reconstruct the level from the position.
 */

void reset_bunchlist()
{	
	#pragma omp for
	for (int i = 0; i < NBunches; i++) { // reconstruct Bunches from Topnodes

		if (D[i].Bunch.Target < 0)
			continue;

		D[i].Bunch.Npart = D[i].Bunch.Cost = D[i].Bunch.Level = 0;
		D[i].Bunch.First_Part = INT_MAX;

		double size = 0.5 * Domain.Size;
		double pos = D[i].TNode.Pos[0];

		while (fmod(pos, size) > 1e-6) {

			D[i].Bunch.Level++;

			size *= 0.5;
		}
	} // for i < NBunches

	return ;
}

/*
 * We split a bunch into 8 sub-bunches/nodes, adding the largest peano key 
 * contained in the bunch. The position of the bunch is set  to a random 
 * particle position during filling. From it we can later construct the 
 * top node center during Tree construction.
 */

static int first_new_bunch = 0;

static int split_bunch(const int parent)
{
	#pragma omp single
	{

	#pragma omp critical // add new nodes at the end
	{

	first_new_bunch = NBunches; 

	NBunches += 8;

	} // omp critical

	for (int i = 0; i < 8; i++) {

		int dest = first_new_bunch + i;

		D[dest].Bunch.Level = D[parent].Bunch.Level + 1;

		int shift = N_SHORT_BITS - 3 * D[dest].Bunch.Level;

		shortKey bitmask = 0x7ULL << shift;
		shortKey keyfragment = ((shortKey) i) << shift;

		D[dest].Bunch.Key = (D[parent].Bunch.Key & ~bitmask) | keyfragment;
		D[dest].Bunch.Npart = 0;
		D[dest].Bunch.First_Part = INT_MAX;
		D[dest].Bunch.Target = -1;
		D[dest].Bunch.Is_To_Be_Split = false;
	}

	} // omp single

	return first_new_bunch;
}

static void remove_empty_bunches()
{
	#pragma omp single
	{

	int old_nBunches = NBunches;

	int i = 0;

	while (i < NBunches) {
	
		if (D[i].Bunch.Npart == 0) {  // remove
	
			NBunches--;

			memmove(&D[i], &D[i+1], (NBunches - i)*sizeof(*D));

			continue;
		} 

		i++;
	}

	} // omp single

	return ;
}

/*
 * update particle distribution over nBunches, starting for first_bunch.
 * Every thread works inside the omp buffer, which are later reduced
 */

static void fill_bunches(const int first_bunch, const int nBunches, 
		 const int first_part, const int nPart)
{
	const int last_part = first_part + nPart;
	const int last_bunch = first_bunch + nBunches;

	struct Bunch_Node *b = Get_Thread_Safe_Buffer(nBunches*sizeof(*b));

	int i = first_bunch;

	for (int j = 0; j < nBunches; j++) {

		b[j].First_Part = INT_MAX;
		b[j].Key = D[i].Bunch.Key;

		i++;
	}

	int j = 0;
	
	#pragma omp for nowait
	for (int ipart = first_part; ipart < last_part; ipart++) {
		
		double px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		double py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		double pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		shortKey pkey = Short_Peano_Key(px, py, pz);

		while (b[j].Key < pkey) // particles are ordered by key
			j++;

		b[j].Npart++;
		b[j].Cost += P[ipart].Cost;
		b[j].First_Part = imin(b[j].First_Part, ipart);
	}

	#pragma omp critical
	{
	
	j = 0;
	
	for (i = first_bunch; i < last_bunch; i++) {

		D[i].Bunch.Npart += b[j].Npart;
		D[i].Bunch.Cost += b[j].Cost;
		D[i].Bunch.First_Part = imin(D[i].Bunch.First_Part, b[j].First_Part);

		j++;
	}

	}// omp critical 

	#pragma omp barrier

	return ;
}

/*
 * This function defines the metric that decides if a bunch has to be refined
 * into eight sub-bunches. It also sets the target processor.
 */

static void distribute()
{
	max_mem_imbal = 0, max_cpu_imbal = 0;

	double mean_npart = Sim.Npart_Total/(Sim.NTask*DOMAIN_NBUNCHES_PER_THREAD);

	Keep_Splitting = false;

	#pragma omp for reduction(max:max_mem_imbal)
	for (int i = 0; i < NBunches; i++ ) {

		double rel_mem_load = (D[i].Bunch.Npart - mean_npart) / mean_npart;

		max_mem_imbal = fmax(max_mem_imbal, rel_mem_load);

		if (rel_mem_load > DOMAIN_SPLIT_MEM_THRES) {

			D[i].Bunch.Is_To_Be_Split = true;
			Keep_Splitting = true;
		}
	}

	return ;
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

static int compare_bunches_by_npart(const void *a, const void *b) 
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Npart > y->Npart) - (x->Npart < y->Npart);
}


static void communicate_particles()
{
	//Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_target);

	return ;
}

static void communicate_bunches()
{
	for (int i = 0; i < NBunches; i++) {

		D[i].Bunch.Target = 0;

		D[i].Bunch.Is_Local = true;
	}
	return ;
}
	

static void print_domain_decomposition (const int max_level)
{
	#pragma omp barrier
	#pragma omp flush(D)

	#pragma omp single
	{
	
	printf(" No | Split | npart  |   sum  | first  | trgt  | lvl ||  PH key\n");
	
	size_t sum = 0;

	for (int i = 0; i < NBunches; i++) {

		sum += D[i].Bunch.Npart;

		printf("%3d |   %d   | %6zu | %6zu | %6d | %5d | %3d || ", 
				i, D[i].Bunch.Is_To_Be_Split, D[i].Bunch.Npart, sum, 
				D[i].Bunch.First_Part, D[i].Bunch.Target, D[i].Bunch.Level);
		
		Print_Int_Bits64(D[i].Bunch.Key);
	}

	} // omp single

	return ;
}

/*
 * Find the global domain origin and the maximum extent. Not much
 * to do for PERIODIC
 */

double max_x = -DBL_MAX, max_y = -DBL_MAX, max_z = -DBL_MAX, 
	   min_x = DBL_MAX, min_y = DBL_MAX, min_z = DBL_MAX;

static void find_global_domain()
{
	Find_Global_Center_Of_Mass(&Domain.Center_Of_Mass[0]);

#ifdef PERIODIC
	
	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;

	Domain.Size = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));

	for (int i = 0; i < 3; i++)
		Domain.Center[i] = Domain.Origin[i] + Domain.Size/2;

#else // ! PERIODIC

	max_x = max_y = max_z = -DBL_MAX;
   	min_x = min_y = min_z = DBL_MAX;

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

	#pragma omp single 
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
#endif // ! PERIODIC

#ifdef DEBUG
	rprintf("\nDomain size %g, \n   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Com    at x = %4g, y = %4g, z = %4g. \n", Domain.Size,
			Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center_Of_Mass[0], Domain.Center_Of_Mass[1],
			Domain.Center_Of_Mass[2]);
#endif // DEBUG

	return ;
}

static double com_x = 0, com_y = 0, com_z = 0, m = 0;

void Find_Global_Center_Of_Mass(double *CoM_out)
{
	com_x = com_y = com_z = m = 0;

	#pragma omp for reduction(+:com_x,com_y,com_z,m) 
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
	
		com_x += P[ipart].Mass * P[ipart].Pos[0];
		com_y += P[ipart].Mass * P[ipart].Pos[1];
		com_z += P[ipart].Mass * P[ipart].Pos[2];
		
		m += P[ipart].Mass;
	}

	double global_com[3] = { com_x, com_y, com_z  };
	double global_m = m;

	#pragma omp single 
	{
		MPI_Allreduce(MPI_IN_PLACE, &global_com, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);
		
		MPI_Allreduce(MPI_IN_PLACE, &global_m, 1, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	} // omp single

	#pragma omp flush

	for (int i = 0; i < 3; i++) 
		CoM_out[i] = global_com[i] / global_m;

	return ;

}
