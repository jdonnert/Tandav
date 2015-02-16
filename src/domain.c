#include "globals.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "peano.h"

static void reconstruct_complete_bunch_list();
static void find_global_domain_extend();
static void fill_bunches(const int, const int, const int, const int);
static void find_mean_cost();
static int remove_empty_bunches();
static void split_bunch(const int, const int);
static void reallocate_topnodes(); // not thread safe
static int check_distribution(const bool);
static void remove_excess_bunches(int *, int *);
static bool distribute();
static void merge_bunch (const int, const int);
static void communicate_particles();
static void communicate_bunches();
static int cost_metric(const int ipart);
static int compare_bunches_by_key(const void *a, const void *b);
static int compare_bunches_by_target(const void *a, const void *b);
static int compare_bunches_by_cost(const void *a, const void *b);
static int compare_load(const void *a, const void *b);
static void print_domain_decomposition (const int);

union Domain_Node_List *D = NULL;

static int Max_NBunches = 0, NTop_Leaves = 0;
static int *Npart = NULL, *Split_Idx = NULL;

static double Top_Node_Alloc_Factor = 0;
static double max_mem_imbal = 0, max_cost_imbal = 0;
static float *Cost = NULL;

static double mean_cost = 0, mean_npart = 0;
/* 
 * Distribute particles in bunches, which are continuous on the Peano curve,
 * but at different level in the tree. The bunches correspond to nodes of 
 * the tree. Bunches also give the top nodes in the tree, with some info 
 * added. 
 * We keep a global list of all the bunches that also contains the 
 * workload and memory footprint of each bunch. An optimal way of 
 * distributing bunches minimises memory and workload imbalance over all
 * Tasks. 
 * To achieve this, we measure mem & cpu cost and refine bunches at the border
 * between tasks: "Split_Idx".  Then we "distribute" them top to bottom across 
 * MPI ranks. This way particle communication is minimised and we avoid the
 * big particle shuffle.
 * Upon reentry we reconstruct the bunchlist to cover the whole domain by 
 * completing the Peano key on every level separately.
 */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

	find_global_domain_extend();

	Sort_Particles_By_Peano_Key();

	reconstruct_complete_bunch_list();

	Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

	fill_bunches(0, NBunches, 0, Task.Npart_Total);

	find_mean_cost();

	int max_level = 0;

	for (;;) {

		#pragma omp single
		remove_empty_bunches();

		Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

		communicate_bunches();

		bool all_done = true;

		#pragma omp single copyprivate(all_done)
		{

		all_done = distribute();

		max_mem_imbal = max_cost_imbal = -DBL_MAX;

		}

		if ((check_distribution(all_done) == 0) && all_done)
			break;

print_domain_decomposition(max_level);

		int old_nBunches = NBunches;

		for (int i = 0; i < old_nBunches; i++ ) {

			if (D[i].Bunch.Modify != 0) { // split into 8

				int first_new_bunch = NBunches;

				#pragma omp single
				if (NBunches + 8 >= Max_NBunches) // make more space !
					reallocate_topnodes();

				split_bunch(i, first_new_bunch);

				fill_bunches(first_new_bunch, 8, D[i].Bunch.First_Part,
						D[i].Bunch.Npart);

				#pragma omp single
				D[i].Bunch.Npart = 0; // mark for deletion
			} // if 
		} // for i
	} // forever

print_domain_decomposition(max_level);
	int nMerged = 0;

	remove_excess_bunches(&nMerged, &max_level);

	rprintf("\nDomain: %d Top Nodes, %d Top Leaves, max level %d  merged %d\n"
			"        Imbalance: Mem %g, Cost %g \n\n", NBunches, NTop_Leaves,
			max_level, nMerged, max_mem_imbal, max_cost_imbal);

//#ifdef DEBUG
	print_domain_decomposition(max_level);
//#endif

	communicate_particles();

	Sig.Tree_Update = true;

	Profile("Domain Decomposition");

	return ;
}

/*
 * Make room for some bunches and build the first node manually.
 */

void Setup_Domain_Decomposition()
{
	Cost = Malloc(Sim.NTask * sizeof(Cost), "Cost");
	Npart = Malloc(Sim.NTask * sizeof(Npart), "Npart");
	Split_Idx = Malloc(Sim.NTask * sizeof(*Split_Idx), "Split_Idx");

	Top_Node_Alloc_Factor = (double) 1024 / Task.Npart_Total;

	reallocate_topnodes();

	memset(D, 0, Max_NBunches * sizeof(*D));

	NBunches = 1;

	D[0].Bunch.Key = 0xFFFFFFFFFFFFFFFF;
	D[0].Bunch.Npart = D[0].Bunch.Level = D[0].Bunch.Target = 0;

	#pragma omp parallel
	{

	find_global_domain_extend();

	rprintf("\nInitial Domain size is %g, \n"
			"   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Center at x = %4g, y = %4g, z = %4g. \n"
			"   CoM    at x = %4g, y = %4g, z = %4g. \n",
			Domain.Size, Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center[0], Domain.Center[1], Domain.Center[2],
			Domain.Center_Of_Mass[0], Domain.Center_Of_Mass[1],
			Domain.Center_Of_Mass[2]);

	} // omp parallel

	return;
}

/*
 * This increases the room for Bunches/Topnodes by 20 %, 
 * so we can stay minimal in memory. Not thread safe ! 
 */

static void reallocate_topnodes()
{
	Top_Node_Alloc_Factor *= 1.2;

	Max_NBunches = Sim.Npart_Total * Top_Node_Alloc_Factor;

	size_t nBytes = Max_NBunches * sizeof(*D);

	printf("Increasing Top Node Memory to %g KB, Max %d Nodes, Factor %4g \n"
			, nBytes/1024.0, Max_NBunches, Top_Node_Alloc_Factor);

	D = Realloc(D, nBytes, "D");

	return ;
}

/*
 * Transform the top nodes back into a bunch list. Add nodes so 
 * the complete domain is covered again. This is equivalent to completing every 
 * triplet from up to 7 (111) starting at the lowest level
 * until the common part of the original keys (top) is reached.
 * Then every triplet of i+1 is incremented with decreasing level until the 
 * bkey triplet is reached.
 *
 *      akey  000.010.010.111.111.111  ->  bkey 010.001.110.001.111.111
 *
 *	000.010.010.111.111.111 > 000.111.111.111.111.111 > 010.000.111.111.111.111
 *	      <- ^                 ^                         ^
 *      <- start 010->111    fill top triplet           start -> 000 -> 001
 *
 * fill every triplet up to 111, fill triplet "top",  fill triplets downwards
 * All new keys are kept in the OmP buffer and later copied back 
 * into "D". At last, we reset properties overwritten by the union in D.
 */

void reconstruct_complete_bunch_list()
{
	const int nOld_Bunches = NBunches;

	if (NBunches == 1)
		goto skip_reconstruction;

	#pragma omp single
	{

	Free(Tree);

	Tree = NULL;

	if (D[NBunches-1].Bunch.Key != 0xFFFFFFFFFFFFFFFF) { // make end

		int i = NBunches++;

		D[i].Bunch.Level = 1;
		D[i].Bunch.Key = 0xFFFFFFFFFFFFFFFF;
		D[i].Bunch.Target = -INT_MAX;
	}

	} // omp single

	struct Bunch_Node *b = Get_Thread_Safe_Buffer(Task.Buffer_Size);

	int nNew = 0;

	#pragma omp for nowait
	for (int i = 0; i < nOld_Bunches-1; i++) { // fill to cover whole domain

		shortKey akey = D[i].Bunch.Key;	// lowest key
		shortKey bkey = D[i+1].Bunch.Key; // highest key

		int top = 1; // highest level where akey != bkey

		uint64_t mask = 0x7ULL << (N_SHORT_BITS-3);

		while ((akey & mask) == (bkey & mask)) {

			top++;
			mask >>= 3;
		}

		for (int j = D[i].Bunch.Level; j > top; j--) { // fill akey upwards

			int shift = N_SHORT_BITS - 3 * j;

			uint64_t atriplet = (akey >> shift) & 0x7;

			uint64_t template = (akey | (0xFFFFFFFFFFFFFFFF >> 3*j))
												& ~(0x7ULL << shift);

			for (uint64_t k = atriplet + 1; k < 8; k++) {

				b[nNew].Key = template | (k << shift);
				b[nNew].Level = j;

				nNew++;
			}
		}

		int shift = N_SHORT_BITS - 3 * top; // fill at level 'top'

		uint64_t atriplet = (akey >> shift) & 0x7;
		uint64_t btriplet = (bkey >> shift) & 0x7;

		uint64_t template = (akey | (0xFFFFFFFFFFFFFFFF >> 3*top))
											 & ~(0x7ULL << shift);

		for (uint64_t k = atriplet + 1; k < btriplet; k++) {

			b[nNew].Level = top;
			b[nNew].Key = template | (k << shift);

			nNew++;
		}

		for (int j = top+1; j <= D[i+1].Bunch.Level; j++) { // fill downwards

			int shift = N_SHORT_BITS - 3 * j;

			uint64_t btriplet = (bkey >> shift) & 0x7;

			uint64_t template = (bkey | (0xFFFFFFFFFFFFFFFF >> 3*j))
														& ~(0x7ULL << shift);
			for (uint64_t k = 0; k < btriplet; k++) {

				b[nNew].Level = j;
				b[nNew].Key = template | (k << shift);

				nNew++;
			} // for k
		} // for j
	} // for i

	int start = 0, end = 0;

	#pragma omp critical
	{

	start = NBunches;
	end = start + nNew;

	while (end >= Max_NBunches)
		reallocate_topnodes();

	NBunches += nNew;

	} // omp critical

	#pragma omp barrier

	int j = 0;

	for (int i = start; i < end; i++) {

		D[i].Bunch.Key = b[j].Key;
		D[i].Bunch.Target = -INT_MAX;
		D[i].Bunch.Level = b[j++].Level;
	}

	skip_reconstruction:;

	#pragma omp for
	for (int i = 0; i < NBunches; i++) { // reset values in all bunches

		D[i].Bunch.Npart = 0;
		D[i].Bunch.Modify = 0;
		D[i].Bunch.Cost = 0;
		D[i].Bunch.First_Part = INT_MAX;

		if (D[i].Bunch.Target >= 0)
			D[i].Bunch.Is_Local = true;

	} // for i < NBunches

	rprintf("Domain: Reconstruction %d -> %d bunches\n", nOld_Bunches,
			NBunches);

	return ;
}

/*
 * Find topnodes to merge, because they are on the same Rank & level 
 * and are complete. We do this until there is nothing left to merge. 
 */

static int nMerged = 0;

static void remove_excess_bunches(int *nMerged_total, int *max_level)
{
	if (NBunches <= Sim.NTask*4)
		return ;

	*nMerged_total = 0;

	nMerged = 1;

	while (nMerged != 0) { // merge recursively

		nMerged = 0;

		#pragma omp single
		//#pragma omp for reduction(+:nMerged)
		for (int i = 0; i < Sim.NTask; i++) {

			int target = -i - 1;
			int first = 1;

			while (D[first].Bunch.Target != target) // find first of target i
				first++;
printf("start i %d target %d first %d \n", i, target, first);

			while (first < NBunches-1) {

				int curr_level = D[first].Bunch.Level;
				int last = 0;

				for (last = first+1; last < NBunches; last++) {

					if (D[last].Bunch.Level != curr_level)
						break;
				}
printf("   currlvl %d first %d last %d \n", curr_level,first,last );

				if (D[last].Bunch.Target != target)
					break;

				if (D[last].Bunch.Level < curr_level
				 && D[first-1].Bunch.Level < curr_level) { // merge
printf("   Merger first %d last %d \n",first, last);
					#pragma omp critical
					merge_bunch(first, last);

					nMerged += last-first;
				}

				first = last;

			} // while target
		} // omp for target

		#pragma omp single
		{

		*max_level = remove_empty_bunches();
		*nMerged_total += nMerged;

		} // omp single

	} // while nMerged

	return ;
}

static void merge_bunch(const int first, const int last)
{
	D[first].Bunch.Level--;
	D[first].Bunch.Modify = 0;
	D[first].Bunch.Key |= 0xFFFFFFFFFFFFFFFF >> (3*D[first].Bunch.Level);

	for (int j = first+1; j < last; j++) { // collapse into first

		D[first].Bunch.Npart += D[j].Bunch.Npart;
		D[first].Bunch.Cost += D[j].Bunch.Cost;

		D[j].Bunch.Npart = 0; // mark for removal
	}

	return ;
}


/*
 * We split a bunch into 8 sub-bunches/nodes, adding the largest peano key 
 * contained in the bunch. The position is later reconstructed from the first 
 * particle contained in the bunch.
 */

static void split_bunch(const int parent, const int first)
{
	#pragma omp single
	NBunches += 8;

	#pragma omp for
	for (int i = 0; i < 8; i++) {

		int dest = first + i;

		D[dest].Bunch.Level = D[parent].Bunch.Level + 1;

		int shift = N_SHORT_BITS - 3 * D[dest].Bunch.Level;

		shortKey bitmask = 0x7ULL << shift;
		shortKey keyfragment = ((shortKey) i) << shift;

		D[dest].Bunch.Key = (D[parent].Bunch.Key & ~bitmask) | keyfragment;
		D[dest].Bunch.Npart = 0;
		D[dest].Bunch.First_Part = INT_MAX;
		D[dest].Bunch.Target = -1;
		D[dest].Bunch.Modify = false;
	}

	return ;
}

static int remove_empty_bunches()
{
	int i = 0;

	int n = NBunches, nLeaves = 0, max_lvl = -1;

	while (i < n) {

		if (D[i].Bunch.Npart == 0) {  // remove

			n--;

			memmove(&D[i], &D[i+1], (n-i) * sizeof(*D)); // fine for n == i 

			continue;
		}

		if (D[i].Bunch.Npart <= 8)
			nLeaves++;

		max_lvl = imax(max_lvl, D[i].Bunch.Level);

		i++;
	}

	NTop_Leaves = nLeaves;
	NBunches = n;

	return max_lvl;
}

/*
 * Update particle distribution over NBunches, starting from first_bunch. 
 * This is performance critical. Every thread works inside its omp buffer, 
 * which are later reduced. The reduction is overlapped with the filling.
 */

static void fill_bunches(const int first_bunch, const int nBunches,
		 const int first_part, const int nPart)
{
	const int last_part = first_part + nPart;
	const int last_bunch = first_bunch + nBunches;

	struct Bunch_Node *b = Get_Thread_Safe_Buffer(nBunches * sizeof(*b));

	int run = first_bunch;

	for (int i = 0; i < nBunches; i++) { // init omp buffer

		b[i].First_Part = INT_MAX;
		b[i].Key = D[run].Bunch.Key;

		run++;
	}

	run = 0;

	#pragma omp for
	for (int i = first_bunch; i < last_bunch; i++) { // reset D

		D[i].Bunch.Npart = D[i].Bunch.Cost = 0;
		D[i].Bunch.First_Part = INT_MAX;
	}

	#pragma omp flush

	#pragma omp for nowait
	for (int ipart = first_part; ipart < last_part; ipart++) { // sort in

		shortKey pkey = Short_Peano_Key(P[ipart].Pos);

		while (b[run].Key < pkey) // particles are ordered by key
			run++;

		b[run].Npart++;
		b[run].Cost += cost_metric(ipart);
		b[run].First_Part = imin(b[run].First_Part, ipart);
	}

	#pragma omp critical
	{

	run = 0;

	for (int i = first_bunch; i < last_bunch; i++) { // reduce

		D[i].Bunch.Npart += b[run].Npart;
		D[i].Bunch.Cost += b[run].Cost;
		D[i].Bunch.First_Part = imin(D[i].Bunch.First_Part,b[run].First_Part);

		run++;
	}

	} // omp critical 

	#pragma omp barrier

	return ;
}

static int cost_metric(const int ipart)
{
	return 1;
}

static void find_mean_cost()
{
	#pragma omp single // do cost and npart mean
	mean_cost = 0;

	#pragma omp for reduction(+:mean_cost)
	for (int i = 0; i < NBunches; i++)
		mean_cost += D[i].Bunch.Cost;

	#pragma omp single
	{

	mean_cost /= Sim.NTask;
	mean_npart = ((double) Sim.Npart_Total) / Sim.NTask;

	}

	return ;
}

/*
 * This function defines the algorithm that decides if a bunch has to be refined
 * into eight sub-bunches. It sets the target processor setting 
 * "target = -rank -1".
 * We refine the bunches adjacent to the processor borders at "Split_Idx".
 */


static int check_distribution(const bool all_done)
{
	if (NBunches == 1 && (Sim.NTask != 1)) { // always split the first

		D[0].Bunch.Modify = 1;

		return 1;
	}

	if (NBunches < Sim.NTask * 4){ // split all

		for (int i = 0; i < NBunches; i++)
			D[i].Bunch.Modify = 2;

		return NBunches;
	}

	#pragma omp for reduction(max:max_mem_imbal,max_cost_imbal) // imbalance
	for (int i = 0; i < Sim.NTask; i++ ) {

		double cost_imbal = fabs(Cost[i] - mean_cost) / mean_cost;

		if (cost_imbal > max_cost_imbal)
			max_cost_imbal = cost_imbal;

		double npart_imbal = fabs(Npart[i] - mean_npart) / mean_npart;

		if (npart_imbal > max_mem_imbal)
			max_mem_imbal = npart_imbal;
		printf("%d %g %g \n",i, cost_imbal, npart_imbal);
	}

	if ((max_mem_imbal < PART_ALLOC_FACTOR-1) // distribution OK ?
	&& (max_cost_imbal < DOMAIN_IMBAL_CEIL))
		return 0;

	int nSplit = 0;

	#pragma omp single copyprivate(nSplit) // put split
	{

	int i = 0;

	for (int task = 0; task < Sim.NTask; task++) {

		if ((Cost[task]/mean_cost > DOMAIN_SPLIT_THRES || ! all_done)
			&& Split_Idx[task] > -1) {

			D[Split_Idx[task]].Bunch.Modify = 1;

			nSplit++;
		}

	} // for task

	} // omp single

	return nSplit;
}

/*
 * Assign tasks to bunches, top to bottom and measure cost.
 */

static bool distribute()
{
	const double max_npart = mean_npart * PART_ALLOC_FACTOR; // per task

	memset(Cost, 0, Sim.NTask * sizeof(*Cost));
	memset(Npart, 0, Sim.NTask * sizeof(*Npart));
	memset(Split_Idx, -1, Sim.NTask * sizeof(*Split_Idx));

	int task = 0;
	int i = 0;
	bool all_done = true;

	for (i = 0; i < NBunches; i++) {

		Cost[task] += D[i].Bunch.Cost;
		Npart[task] += D[i].Bunch.Npart;

		D[i].Bunch.Target = -task - 1;

		if ((Cost[task] > mean_cost) || (Npart[task] > max_npart)) {

			if (D[i].Bunch.Cost > D[i+1].Bunch.Cost)
				Split_Idx[task] = i;
			else
				Split_Idx[task] = i+1;

			task++;
		}

		if (task == Sim.NTask) {

			all_done = false;

			Split_Idx[Sim.NTask-1] = i;

			break;
		}
	} // for i


	return all_done;
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

static int compare_bunches_by_cost(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Cost > y->Cost) - (x->Cost < y->Cost);
}

static void communicate_particles()
{
	//Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_target);

	return ;
}

/*
 * Reduce the Bunch list over all MPI ranks 
 */

static void communicate_bunches()
{
	#pragma omp for
	for (int i = 0; i < NBunches; i++) {

		//D[i].Bunch.Target = 0;

		D[i].Bunch.Is_Local = true;
	}


	// MPI_Ibcast();

	return ;
}



/*
 * Find the global domain origin and the maximum extent. We center the domain 
 * on the center of mass to make the decomposition effectively Lagrangian if
 * the simulation is not PERIODIC. For PERIODIC simulations there is little to
 * do.
 */

double max_distance = 0;

static void find_global_domain_extend()
{
	Find_Global_Center_Of_Mass(&Domain.Center_Of_Mass[0]);

#ifdef PERIODIC

	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;

	Domain.Size = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));

	for (int i = 0; i < 3; i++)
		Domain.Center[i] = Domain.Origin[i] + 0.5 * Domain.Size;

#else // ! PERIODIC

	#pragma omp single
	max_distance = 0;

	#pragma omp for reduction(max:max_distance)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		for (int i = 0; i < 3; i++) {

			if (P[ipart].Pos[i] > max_distance)
				max_distance = P[ipart].Pos[i];

			if (-1*P[ipart].Pos[i] > max_distance)
				max_distance = -1*P[ipart].Pos[i];

		} // for i
	} // for ipart

	#pragma omp single
	{

	MPI_Allreduce(MPI_IN_PLACE, &max_distance, 1, MPI_DOUBLE, MPI_MAX,
		MPI_COMM_WORLD);

	Domain.Size = 2.0 * max_distance;

	for (int i = 0; i < 3; i++) {

		Domain.Origin[i] = - 0.5 * Domain.Size;
		Domain.Center[i] = Domain.Origin[i] + 0.5 * Domain.Size ;
	}

	} // omp single (Domain)

	#pragma omp flush

#endif // ! PERIODIC

#ifdef DEBUG
	rprintf("\nDomain size is %g, \n"
			"   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Center at x = %4g, y = %4g, z = %4g. \n"
			"   CoM    at x = %4g, y = %4g, z = %4g. \n",
			Domain.Size, Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center[0], Domain.Center[1], Domain.Center[2],
			Domain.Center_Of_Mass[0], Domain.Center_Of_Mass[1],
			Domain.Center_Of_Mass[2]);
#endif

	return ;
}

static double com_x = 0, com_y = 0, com_z = 0, m = 0;

void Find_Global_Center_Of_Mass(double *CoM_out)
{
	#pragma omp single
	com_x = com_y = com_z = m = 0;

	#pragma omp barrier

	#pragma omp for reduction(+:com_x,com_y,com_z,m)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		com_x += P[ipart].Mass * P[ipart].Pos[0];
		com_y += P[ipart].Mass * P[ipart].Pos[1];
		com_z += P[ipart].Mass * P[ipart].Pos[2];

		m += P[ipart].Mass;
	}

	#pragma omp single
	{

	double global_com[3] = { com_x, com_y, com_z  };
	double global_m = m;

	MPI_Allreduce(MPI_IN_PLACE, &global_com, 3, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, &global_m, 1, MPI_DOUBLE, MPI_MIN,
			MPI_COMM_WORLD);

	for (int i = 0; i < 3; i++)
		CoM_out[i] = global_com[i] / global_m;

	} // omp single

	return ;

}

static void print_domain_decomposition (const int max_level)
{

	#pragma omp flush

	#pragma omp master
	{

	rprintf(" No | Split | npart  |   sum  | first  | trgt  | lvl |  Cost  | "
			"CumCost | Max PH key,   Max_level %d \n", max_level);

	size_t sum = 0;
	double csum = 0;
	int last_target = -1;

	for (int i = 0; i < NBunches; i++) {

		sum += D[i].Bunch.Npart;

		if (last_target != D[i].Bunch.Target){

			last_target = D[i].Bunch.Target;
			csum = 0;
		}

		csum += D[i].Bunch.Cost;

		printf("%3d |   %d   | %6zu | %6zu | %6d | %5d | %3d | %6g | %6g || ",
				i, D[i].Bunch.Modify, D[i].Bunch.Npart, sum,
				D[i].Bunch.First_Part, D[i].Bunch.Target, D[i].Bunch.Level,
				D[i].Bunch.Cost, csum);

		Print_Int_Bits64(D[i].Bunch.Key);
	}

	Assert(sum == Sim.Npart_Total, "More or less particles in D than in Sim");

	} // omp master

	#pragma omp barrier

	return ;
}
