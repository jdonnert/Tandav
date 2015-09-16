#include "globals.h"
#include "Gravity/gravity.h"
#include "domain.h"
#include "peano.h"

static void set_global_domain();
static void fill_bunches(const int, const int, const int, const int);
static void find_mean_cost();
static int remove_empty_bunches();
static void split_bunch(const int, const int);
static void reallocate_topnodes(); // not thread safe
static int check_distribution();
static void remove_excess_bunches();
static void reset_bunchlist();
static void distribute();
static void merge_bunch (const int, const int);
static void communicate_particles();
static void communicate_bunches();
static int cost_metric(const int ipart);
static int compare_bunches_by_key(const void *a, const void *b);
static int compare_bunches_by_target(const void *a, const void *b);
static int compare_bunches_by_cost(const void *a, const void *b);
static void print_domain_decomposition (const int);
static void find_domain_center(double *Center_out);
static double find_largest_particle_distance();

union Domain_Node_List *D = NULL;

static int NTarget = 0;
static int Max_NBunches = 0, NTop_Leaves = 0;
static int NMerged = 0, Max_Level = 0;
static int *Npart = NULL, *Split_Idx = NULL;

static double max_mem_imbal = 0, max_cost_imbal = 0;
static float *Cost = NULL;

static double Mean_Cost = 0, Mean_Npart = 0;

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
 * Upon reentry we start off with one top node only, as this is a log(n)
 * algorithm.
 */

void Domain_Decomposition()
{
	Profile("Domain Decomposition");

printf("DOMAIN %d \n", Task.Thread_ID); fflush(stdout);
	set_global_domain();

printf("SORT %d \n", Task.Thread_ID);fflush(stdout);
	Sort_Particles_By_Peano_Key();

printf("RESET %d \n", Task.Thread_ID);fflush(stdout);
	reset_bunchlist();

printf("FILL %d \n", Task.Thread_ID);fflush(stdout);
	fill_bunches(0, NBunches, 0, Task.Npart_Total);

	find_mean_cost();

	for (;;) {

		#pragma omp single
		remove_empty_bunches();

		Qsort(Sim.NThreads, D, NBunches, sizeof(*D), &compare_bunches_by_key);

		communicate_bunches();

		#pragma omp single
		{

		distribute();

		max_mem_imbal = max_cost_imbal = -DBL_MAX;

		}

		if ((check_distribution() == 0))
			break;

		int old_nBunches = NBunches;

		for (int i = 0; i < old_nBunches; i++ ) {

			if (D[i].Bunch.Modify != 0) { // split into 8

				int first_new_bunch = NBunches;

				#pragma omp single
				if (NBunches + 8 >= Max_NBunches) // make more space !
					reallocate_topnodes();
printf("SPLIT %d \n", Task.Thread_ID);
				split_bunch(i, first_new_bunch);

printf("FILL %d \n", Task.Thread_ID);
				fill_bunches(first_new_bunch, 8, D[i].Bunch.First_Part,
						D[i].Bunch.Npart);

				#pragma omp single
				D[i].Bunch.Npart = 0; // mark for deletion
			} // if 
printf("OUT %d \n", Task.Thread_ID);
		} // for i
	} // forever

	remove_excess_bunches();

	rprintf("\nDomain: %d Top Nodes, %d Top Leaves, max level %d  merged %d\n"
			"        Max Imbalance: Mem %g, Cost %g \n\n", NBunches,
			NTop_Leaves, Max_Level, NMerged, max_mem_imbal, max_cost_imbal);

#ifdef DEBUG
	print_domain_decomposition(Max_Level);
#endif

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
	if (Sim.NRank < 16) // decompose on threads as welL
		NTarget = Sim.NTask;
	else
		NTarget = Sim.NRank;

	Cost = Malloc(Sim.NTask * sizeof(Cost), "Domain Cost");
	Npart = Malloc(Sim.NTask * sizeof(Npart), "Domain Npart");
	Split_Idx = Malloc(Sim.NTask * sizeof(*Split_Idx), "Domain Split_Idx");

	int min_level = log(Sim.NTask)/log(8) + 1;

	Max_NBunches = pow(8, min_level);

	reallocate_topnodes();

	reset_bunchlist();

	#pragma omp parallel
	set_global_domain();

	rprintf("\nDomain size is %g, \n"
		"   Origin at x = %4g, y = %4g, z = %4g, \n"
		"   Center at x = %4g, y = %4g, z = %4g. \n"
		"   CoM    at x = %4g, y = %4g, z = %4g. \n",
		Domain.Size, Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
		Domain.Center[0], Domain.Center[1], Domain.Center[2],
		Domain.Center_Of_Mass[0], Domain.Center_Of_Mass[1],
		Domain.Center_Of_Mass[2]); fflush(stdout);

	return;
}

/*
 * This increases the room for Bunches/Topnodes by 20 %, 
 * so we can stay minimal in memory. Not thread safe ! 
 */

static void reallocate_topnodes()
{
	Max_NBunches *= 1.2;

	size_t nBytes = Max_NBunches * sizeof(*D);

	double alloc_factor = Max_NBunches * sizeof(*D) /
						((double) Sim.Npart_Total* sizeof(*P));

	printf("Increasing Top Node Memory by 20%% to %g KB, Max %d Nodes, "
		   "Factor %4g \n",
		   nBytes/1024.0, Max_NBunches, alloc_factor);

	D = Realloc(D, nBytes, "D");

	return ;
}


/*
 * We set a vanilla bunch list with a minimum level to allow immediate omp
 * parallel fill.
 */

static void reset_bunchlist()
{
	#pragma omp single
	memset(&D[0], 0, sizeof(*D) * Max_NBunches);

	int level = log(Sim.NTask)/log(8) + 1;
	NBunches = pow(8,level);

	int shift = 3 * level;
	shortKey base = 0xFFFFFFFFFFFFFFFF >> shift;

	#pragma omp for
	for (shortKey i = 0; i < NBunches; i++) {

		D[i].Bunch.Key = (i << (64 - shift)) | base;
		D[i].Bunch.Level = level;
		D[i].Bunch.Npart = 0;
		D[i].Bunch.First_Part = INT_MAX;
		D[i].Bunch.Target = -1;
		D[i].Bunch.Modify = 0;
	}

	return ;
}

/*
 * Find topnodes to merge, because they are on the same Rank & level 
 * and are complete. We do this until there is nothing left to merge. 
 */

static void remove_excess_bunches()
{
	if (NBunches <= NTarget*2)
		return ;

	#pragma omp single
	NMerged = Max_Level = 0;

	#pragma omp for reduction(+:NMerged) reduction(max:Max_Level)
	for (int i = 0; i < Sim.NTask - 1; i++) {

		int target = -i - 1;

		int split = Split_Idx[i];
		int split_level = D[split].Bunch.Level;

		int last = split - 1;
		int first = last;

		while (D[first].Bunch.Level > split_level)
			first--;

		first++;

		if (first >= last)
			continue;

		if (D[first].Bunch.Target != target)
			continue;

		D[first].Bunch.Level = split_level;
		D[first].Bunch.Modify = 2;
		D[first].Bunch.Key |= 0xFFFFFFFFFFFFFFFF >> (3 * split_level);

		for (int j = first+1; j <= last; j++) { // collapse into first

			D[first].Bunch.Npart += D[j].Bunch.Npart;
			D[first].Bunch.Cost += D[j].Bunch.Cost;

			D[j].Bunch.Npart = 0; // mark for removal
			D[j].Bunch.Modify = 2;
		}

		NMerged += last-first;

	} // omp for target

	#pragma omp single
	Max_Level = remove_empty_bunches();

	return ;
}

static void merge_bunch(const int first, const int last)
{
	D[first].Bunch.Level--;
	D[first].Bunch.Modify = 2;
	D[first].Bunch.Key |= 0xFFFFFFFFFFFFFFFF >> (3*D[first].Bunch.Level);

	for (int j = first+1; j <= last; j++) { // collapse into first

		D[first].Bunch.Npart += D[j].Bunch.Npart;
		D[first].Bunch.Cost += D[j].Bunch.Cost;

		D[j].Bunch.Npart = 0; // mark for removal
		D[j].Bunch.Modify = 2;
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
		D[dest].Bunch.Modify = 0;
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

	struct Bunch_Node *buf = Get_Thread_Safe_Buffer(nBunches * sizeof(*buf));

	int run = first_bunch;

	for (int i = 0; i < nBunches; i++) { // init omp buffer

		buf[i].First_Part = INT_MAX;
		buf[i].Key = D[run].Bunch.Key;

		run++;
	}

	run = 0;

	#pragma omp for
	for (int i = first_bunch; i < last_bunch; i++) { // reset D

		D[i].Bunch.Npart = D[i].Bunch.Cost = 0;
		D[i].Bunch.First_Part = INT_MAX;
	}

	#pragma omp for nowait
	for (int ipart = first_part; ipart < last_part; ipart++) { // sort in

		shortKey pkey = Short_Peano_Key(P[ipart].Pos);

		while (buf[run].Key < pkey) // particles are ordered by key
			run++;

		buf[run].Npart++;
		buf[run].Cost += cost_metric(ipart);
		buf[run].First_Part = imin(buf[run].First_Part, ipart);
	}

	#pragma omp critical
	{

	run = 0;

	for (int i = first_bunch; i < last_bunch; i++) { // reduce

		D[i].Bunch.Npart += buf[run].Npart;
		D[i].Bunch.Cost += buf[run].Cost;
		D[i].Bunch.First_Part = imin(D[i].Bunch.First_Part,
									 buf[run].First_Part);
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
	Mean_Cost = 0;

	#pragma omp for reduction(+:Mean_Cost)
	for (int i = 0; i < NBunches; i++)
		Mean_Cost += D[i].Bunch.Cost;

	#pragma omp single
	{

	Mean_Cost /= NTarget;
	Mean_Npart = ((double) Sim.Npart_Total) / NTarget;

	}

	return ;
}

/*
 * This function defines the algorithm that decides if a bunch has to be 
 * refined into eight sub-bunches. It sets the target processor setting 
 * "target = -rank -1".
 * We refine the bunches adjacent to the processor borders at "Split_Idx".
 */


static int check_distribution()
{
	if (NBunches == 1) { // always split the first

		D[0].Bunch.Modify = 1;

		return 1;
	}

	#pragma omp for reduction(max:max_mem_imbal,max_cost_imbal) // imbalance
	for (int i = 0; i < Sim.NTask; i++ ) {

		double cost_imbal = fabs(Cost[i] - Mean_Cost) / Mean_Cost;

		if (cost_imbal > max_cost_imbal)
			max_cost_imbal = cost_imbal;

		double npart_imbal = fabs(Npart[i] - Mean_Npart) / Mean_Npart;

		if (npart_imbal > max_mem_imbal)
			max_mem_imbal = npart_imbal;
	}

	if ((max_mem_imbal < PART_ALLOC_FACTOR-1) // distribution OK ?
	&& (max_cost_imbal < DOMAIN_IMBAL_CEIL)
	&& (NBunches > Sim.NThreads) )
		return 0;

	int nSplit = 0;

	#pragma omp single copyprivate(nSplit) // put split mark
	{

	for (int task = 0; task < Sim.NTask; task++) {

		if (fabs(Cost[task]-Mean_Cost)/Mean_Cost < DOMAIN_IMBAL_CEIL)
			continue;

		if (Split_Idx[task] < 0) // no idx set
			continue;

		int i = Split_Idx[task];

		if (D[i].Bunch.Level == N_SHORT_TRIPLETS) // maximum refinement
			continue;

		D[i].Bunch.Modify = 1;

		nSplit++;
	} // for task

	} // omp single

	return nSplit;
}

/*
 * Assign tasks to bunches, top to bottom and measure cost.
 * We accept a bunch if it brings us closer to the cost mean and mark the 
 * bunch hitting the mean as to split. We adjust for the imbalance in counting
 */

static void distribute()
{
	double max_npart = Mean_Npart * PART_ALLOC_FACTOR;

	memset(Cost, 0, Sim.NTask * sizeof(*Cost));
	memset(Npart, 0, Sim.NTask * sizeof(*Npart));
	memset(Split_Idx, -1, Sim.NTask * sizeof(*Split_Idx));

	int task = 0;

	double dCost = 0, dNpart = 0;

	for (int i = 0; i < NBunches; i++) {

		double cur_cost = Cost[task] + dCost;
		double new_cost = cur_cost + D[i].Bunch.Cost;
		double new_npart = Npart[task] + dNpart + D[i].Bunch.Npart;

		double dmean_new = fabs(new_cost - Mean_Cost);
		double dmean_cur = fabs(cur_cost - Mean_Cost);

		if ( (new_npart >= max_npart) || dmean_new > dmean_cur ) { // new rank

			if (Cost[task] > Mean_Cost)
				Split_Idx[task] = i-1;
			else
				Split_Idx[task] = i;

			task++;

			int last_task = task - 1;

			task = task % NTarget;
			last_task = last_task % NTarget;

			dCost = Cost[last_task] + dCost - Mean_Cost; // carry the delta
			dNpart = Npart[last_task] + dNpart - Mean_Npart;
		}

		Cost[task] += D[i].Bunch.Cost;
		Npart[task] += D[i].Bunch.Npart;

		D[i].Bunch.Target = -task - 1;
	} // for i

	return ;
}

static int compare_bunches_by_key(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Key > y->Key) - (x->Key < y->Key);
}

static int compare_bunches_by_cost(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Cost > y->Cost) - (x->Cost < y->Cost);
}

static void communicate_particles()
{

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
 * the simulation is not PERIODIC. Actually the median would be much better !
 * For PERIODIC simulations there is little to do.
 */


static void set_global_domain()
{
#ifdef PERIODIC

	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;

	Domain.Size = fmax(Sim.Boxsize[0], fmax(Sim.Boxsize[1], Sim.Boxsize[2]));

#else

	Domain.Size = find_largest_particle_distance();

	find_domain_center(Domain.Center);

#endif // ! PERIODIC

	#pragma omp single
	for (int i = 0; i < 3; i++)
		Domain.Origin[i] = Domain.Center[i] - 0.5 * Domain.Size;

#ifdef DEBUG
	rprintf("\nDomain size is %g, \n"
			"   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Center at x = %4g, y = %4g, z = %4g. \n"
			"   CoM    at x = %4g, y = %4g, z = %4g. \n",
			Domain.Size, Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center[0], Domain.Center[1], Domain.Center[2],
			Sim.Center_Of_Mass[0], Sim.Center_Of_Mass[1],
			Sim.Center_Of_Mass[2]);
#endif

	return ;
}

/*
 * Domain Center is not the center of mass but the median of mass, which is
 * less sensitive to outliers.
 */

static Float *x = NULL;

static void find_domain_center(double Center_out[3])
{
	Float center[3] = { 0 };

	#pragma omp single
	x = Malloc(Task.Npart_Total * sizeof(*x), "x");

	for (int i = 0; i < 3; i++) {

		#pragma omp for
		for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
			x[ipart] = P[ipart].Mass * P[ipart].Pos[i];

		center[i] = Median(Task.Npart_Total, x) / Sim.Total_Mass;

		#pragma omp barrier
	}

	#pragma omp single
	Free(x);

	#pragma omp master
	{

	Float sendbuf[Sim.NRank], recvbuf[Sim.NRank];

	for (int i = 0; i < 3; i++) {

		sendbuf[Task.Rank] = center[i];

		MPI_Gather(recvbuf, 1, MPI_MYFLOAT, sendbuf, 1, MPI_MYFLOAT,
				   Sim.Master, MPI_COMM_WORLD);

		center[i] = Median(Sim.NRank, recvbuf);
	}

	MPI_Scatter(center, 3, MPI_DOUBLE, Center_out, 3, MPI_DOUBLE,
				Sim.Master, MPI_COMM_WORLD);

	} // master

	#pragma omp barrier

	return ;
}

static double Max_Distance = 0;

static double find_largest_particle_distance()
{
	#pragma omp single
	Max_Distance = 0;

	#pragma omp for reduction(max:Max_Distance)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {

		for (int i = 0; i < 3; i++) {

			if (P[ipart].Pos[i] > Max_Distance)
				Max_Distance = P[ipart].Pos[i] - Domain.Center_Of_Mass[i];

			if (-1*P[ipart].Pos[i] > Max_Distance)
				Max_Distance = -1*P[ipart].Pos[i] - Domain.Center_Of_Mass[i];

		} // for i
	} // for ipart

	#pragma omp master
	MPI_Allreduce(MPI_IN_PLACE, &Max_Distance, 1, MPI_DOUBLE, MPI_MAX,
		MPI_COMM_WORLD);


	return 2.001 * Max_Distance; // 2.001 helps with cancellation
}

static void print_domain_decomposition (const int max_level)
{
	#pragma omp master
	{

	printf(" No | Split | npart  |   sum  | first  | trgt  | lvl |  Cost  | "
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
