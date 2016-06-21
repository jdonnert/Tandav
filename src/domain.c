#include "domain.h"

#define MIN_LEVEL 2 // decompose at least 8^MIN_LEVEL domains downward

static void set_computational_domain();
static void find_domain_center(double Center_Out[3]);
static void find_largest_particle_distance(double *);
static void reset_bunchlist();
static void fill_new_bunches(const int, const int, const int, const int);
static void find_mean_cost();
static void make_new_bunchlist();
static int remove_empty_bunches();
static void split_bunch(const int, const int);
static void reallocate_topnodes(); // not thread safe
//static void remove_excess_bunches();
static int find_min_level();
static void transform_bunches_into_top_nodes();
static void distribute();
static void find_global_imbalances();
static void mark_bunches_to_split();
static unsigned int cost_metric(const int ipart);
static void zero_particle_cost();
static int compare_bunches_by_key(const void *a, const void *b);
static int compare_bunches_by_cost(const void *a, const void *b);
static void communicate_particles();
static void communicate_bunches();
static void communicate_top_nodes();

#ifdef DEBUG_DOMAIN
static void print_domain_decomposition(const int);
#else
static inline void print_domain_decomposition(const int a) {};
#endif // DEBUG_DOMAIN


union Domain_Node_List * restrict D = NULL;

static int NTarget = 0;
static int Max_NBunches = 0, NTop_Leaves = 0;
static int NMerged = 0, Max_Level = 0;

static int * restrict Npart = NULL, * restrict Split_Idx = NULL;
static float * restrict Cost = NULL;

static double Max_Mem_Imbal = -DBL_MAX, Max_Cost_Imbal = -DBL_MAX;
static double Mean_Cost = 0, Mean_Npart = 0;

static int NBunches = 0;

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

	set_computational_domain();

	Sort_Particles_By_Peano_Key();

	reset_bunchlist();

	fill_new_bunches(0, NBunches, 0, Task.Npart_Total);
	
	print_domain_decomposition(Max_Level); // DEBUG_DOMAIN

	find_mean_cost();

	int cnt = 0;
	
	for (;;) {

		#pragma omp single
		Max_Level = remove_empty_bunches();

		Qsort(D, NBunches, sizeof(*D), &compare_bunches_by_key);

		communicate_bunches();

		distribute();

		find_global_imbalances();
		
		if (Max_Cost_Imbal < DOMAIN_IMBAL_CEIL)
			if (Max_Mem_Imbal < Param.Part_Alloc_Factor-1) // distrib OK ?
				break;

		if (cnt++ > N_SHORT_TRIPLETS - MIN_LEVEL) {
	
			#pragma omp master
			Warn(true, "Domain Decomposition not optimal !");
			
			break;
		}
		
		#pragma omp single
		mark_bunches_to_split();

		make_new_bunchlist();
		
	} // forever

	//remove_excess_bunches();

	print_domain_decomposition(Max_Level); // DEBUG_DOMAIN

	rprintf("\nDomain: After %d iterations ...\n"
			"        %d Top Nodes, %d Top Leaves, max level %d  merged %d\n"
			"        Max Imbalance: Mem %g, Cost %g \n"
			"        Mean Cost %g, Mean Npart %g \n\n", cnt, 
			NBunches, NTop_Leaves, Max_Level, NMerged, Max_Mem_Imbal, 
			Max_Cost_Imbal, Mean_Cost, Mean_Npart);

	communicate_particles();

	transform_bunches_into_top_nodes();

	communicate_top_nodes();

	Sort_Particles_By_Peano_Key();

	Reverse_Peano_Keys();

	Sig.Tree_Update = true;

	Profile("Domain Decomposition");

	return ;
}

/*
 * Make room for some bunches 
 */

void Setup_Domain_Decomposition()
{
	if (NTask < 65) // decompose on threads as welL
		NTarget = NTask;
	else
		NTarget = NRank;

	Cost = Malloc(NTarget * sizeof(Cost), "Domain Cost");
	Npart = Malloc(NTarget * sizeof(Npart), "Domain Npart");
	Split_Idx = Malloc(NTarget * sizeof(*Split_Idx), "Domain Split_Idx");

	int min_level = find_min_level();

	Max_NBunches = pow(8, min_level);

	reallocate_topnodes();

	#pragma omp parallel
	{

	reset_bunchlist();

	set_computational_domain();

	Sort_Particles_By_Peano_Key();
	

	} // omp parallel

	rprintf("\nDomain size is %g, \n"
			"   Origin at x = %4g, y = %4g, z = %4g, \n"
			"   Center at x = %4g, y = %4g, z = %4g. \n"
			"   CoM    at x = %4g, y = %4g, z = %4g. \n", Domain.Size, 
			Domain.Origin[0], Domain.Origin[1], Domain.Origin[2],
			Domain.Center[0], Domain.Center[1], Domain.Center[2],
			Prop.Center_Of_Mass[0], Prop.Center_Of_Mass[1],
			Prop.Center_Of_Mass[2]);

	return;
}

void Finish_Domain_Decomposition()
{
	Free(Cost);
	Free(Npart);
	Free(Split_Idx);
	
	return ;
}

static void communicate_top_nodes()
{
/*	MPI_Request *request = NULL;
		
		float *target = &D[i].TNode.Pos[0];
		int nBytes = (&D[i].TNode.Dp[2] - target) + sizeof(float);

		MPI_Ibcast(target, nBytes, MPI_BYTE, src, MPI_COMM_WORLD, request); */

	return ;
}



/*
 * Set the "TNode" part of the "D"omain unions. From here onwards the members 
 * are to be understood as a Top Node, not a bunch. Also returns the first
 * particle in "ipart" and the "level" 
 */

static void transform_bunches_into_top_nodes()
{
	#pragma omp single nowait
	NTop_Nodes = NBunches;

	#pragma omp for
	for (int i = 0; i < NBunches; i++) {

		uint64_t npart = D[i].Bunch.Npart;

		Assert(npart < INT_MAX, "Npart %zu in Bunch %d > INT_MAX %d",
		   npart, i, INT_MAX);
		
		D[i].TNode.Npart = npart;

		int ipart = D[i].TNode.First_Part;

		double px = P.Pos[0][ipart] - Domain.Origin[0]; 	
		double py = P.Pos[1][ipart] - Domain.Origin[1];
		double pz = P.Pos[2][ipart] - Domain.Origin[2];

		double size = Domain.Size / (1ULL << D[i].TNode.Level);

		D[i].TNode.Npart = npart;
		D[i].TNode.Pos[0] = (floor(px/size) + 0.5) * size + Domain.Origin[0];
		D[i].TNode.Pos[1] = (floor(py/size) + 0.5) * size + Domain.Origin[1];
		D[i].TNode.Pos[2] = (floor(pz/size) + 0.5) * size + Domain.Origin[2];
	}

	return ;
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
						((double) Sim.Npart_Total* sizeof_P);

	rprintf("Increasing Top Node Memory by 20%% to %g KB, Max %d Nodes, "
		   "Factor %4g \n",
		   nBytes/1024.0, Max_NBunches, alloc_factor);

	D = Realloc(D, nBytes, "D");

	return ;
}


/*
 * Clear bunchlist. Also deallocates the Tree
 */

static void reset_bunchlist()
{
	
	#pragma omp barrier
	
	memset(&D[0], 0, sizeof(*D) * Max_NBunches);
	
	int level = find_min_level();
		
	#pragma omp single
	NBunches = pow(8, level);

	const int shift = 3 * level;
	shortKey basekey = 0xFFFFFFFFFFFFFFFF >> shift;
	
	#pragma omp for 
	for (int i = 0; i < NBunches; i++) {

		shortKey key = i;

		D[i].Bunch.Key = (key << (64 - shift)) | basekey;
		D[i].Bunch.Level = level;
		D[i].Bunch.Npart = 0;
		D[i].Bunch.First_Part = INT_MAX;
		D[i].Bunch.Target = -1;
		D[i].Bunch.Modify = 0;
	}

	return ;
}

static int find_min_level()
{
	return fmax(MIN_LEVEL, log(NTarget)/log(8) + 1);

}

/*
 * Find bunches to merge, because they are on the same Rank & level 
 * and are complete. We do this until there is nothing left to merge. 
 */

/*static void remove_excess_bunches()
{
	if (NBunches <= NTarget*2)
		return ;

	#pragma omp single
	NMerged = Max_Level = 0;

	#pragma omp for reduction(+:NMerged) reduction(max:Max_Level)
	for (int i = 0; i < NTarget - 1; i++) {

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
*/ 

static void make_new_bunchlist()
{
	int old_nBunches = NBunches;

	for (int i = 0; i < old_nBunches; i++ ) {

		if (D[i].Bunch.Modify == 0)
			continue;

		int first_new_bunch = NBunches; // split this one into 8

		#pragma omp barrier

		#pragma omp single
		if (NBunches + 8 >= Max_NBunches) // make more space !
			reallocate_topnodes();

		split_bunch(i, first_new_bunch);

		fill_new_bunches(first_new_bunch, 8, D[i].Bunch.First_Part,
				D[i].Bunch.Npart);

		#pragma omp single
		D[i].Bunch.Npart = 0; // mark for deletion

	} // for i

	return ; 
}

/*
 * We split a bunch into 8 sub-bunches/nodes, adding the largest peano key 
 * contained in the bunch. The position is later reconstructed from the first 
 * particle contained in the bunch, when we transform into top nodes.
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

static void fill_new_bunches(const int first_bunch, const int nBunches,
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

		shortKey pkey = (P.Key[ipart] >> DELTA_PEANO_BITS);

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

static unsigned int cost_metric(const int ipart)
{
	return fmax(1, P.Cost[ipart]);
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

	} // omp single

	return ;
}

static void zero_particle_cost()
{
	#pragma omp for 
	for (int i = 0; i < NBunches; i++)
		D[i].Bunch.Cost = 0;

	return ;
}

/*
 * Assign tasks to bunches, top to bottom and measure cost.
 * We use the standard Gadget way of distributing, which is: order the bunches
 * by cost and then assign CPUs.
 */

static void distribute()
{
	#pragma omp single 
	{
	
	memset(Cost, 0, NTarget * sizeof(*Cost));
	memset(Npart, 0, NTarget * sizeof(*Npart));
	memset(Split_Idx, -1, NTarget * sizeof(*Split_Idx));
	
	}
		
	Qsort(&D[0], NBunches, sizeof(*D), &compare_bunches_by_cost);
	
	#pragma omp single
	for (int i = 0; i < NBunches; i++) {

		int task = 0; 
	
		double max_delta_mean = 0;

		for (int j = 0; j < NTarget; j++) { // find a task for bunch i
		
			double delta_mean = Mean_Cost - Cost[j];

			if (delta_mean > max_delta_mean) {
			
				task = j;
				max_delta_mean = delta_mean;
			}
		}
	
		#pragma omp atomic
		Cost[task] += D[i].Bunch.Cost;
		#pragma omp atomic
		Npart[task] += D[i].Bunch.Npart;

		D[i].Bunch.Target = -task - 1;

	} // for i

	Qsort(&D[0], NBunches, sizeof(*D), &compare_bunches_by_key);

	for (int i = 0; i < NTarget; i++) 
			Split_Idx[i] = i;

	for (int i = 0; i < NBunches; i++) {
	
		int task = -1 * (D[i].Bunch.Target + 1);
		int idx = Split_Idx[task];
		
		if (D[i].Bunch.Cost > D[idx].Bunch.Cost)
			Split_Idx[task] = i;
	}

	return ;
}

static void find_global_imbalances()
{
	#pragma omp single
	Max_Mem_Imbal = Max_Cost_Imbal = -DBL_MAX;

	#pragma omp for reduction(max:Max_Mem_Imbal,Max_Cost_Imbal) // imbalance
	for (int i = 0; i < NTarget; i++ ) {

		double cost_imbal = fabs(Cost[i] - Mean_Cost) / Mean_Cost;

		if (cost_imbal > Max_Cost_Imbal)
			Max_Cost_Imbal = cost_imbal;

		double npart_imbal = fabs(Npart[i] - Mean_Npart) / Mean_Npart;

		if (npart_imbal > Max_Mem_Imbal)
			Max_Mem_Imbal = npart_imbal;
	}

	return ;
}

/*
 * This function decides if a bunch can be refined into eight sub-bunches. 
 */

static void mark_bunches_to_split()
{
	for (int i = 0; i < NTarget; i++) {

		//if (fabs(Cost[i]-Mean_Cost)/Mean_Cost < DOMAIN_IMBAL_CEIL)
		//	continue;

		if (Split_Idx[i] < 0) // no idx set
			continue;

		int j = Split_Idx[i];

		if (D[j].Bunch.Level == N_SHORT_TRIPLETS) // maximum refinement
			continue;

		D[j].Bunch.Modify = 1;
	
	} // for i

	return ;
}

static int compare_bunches_by_key(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Key < y->Key) - (x->Key > y->Key);
}

static int compare_bunches_by_cost(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Cost < y->Cost) - (x->Cost > y->Cost);
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

	// MPI_Allreducev();
	// MPI_Ibcast();

	return ;
}

static void communicate_particles()
{
	zero_particle_cost();

	return ;
}

/*
 * Find the global domain origin and the maximum extent. We center the domain 
 * on the center of mass to make the decomposition effectively Lagrangian if
 * the simulation is not PERIODIC. Actually the median would be much better !
 * For PERIODIC simulations there is little to do.
 */

#ifdef PERIODIC

static void set_computational_domain()
{
	#pragma omp single
	{
	
	Domain.Size = Sim.Boxsize[0];

	Domain.Center[0] = Domain.Center[1] = Domain.Center[2] = 0.5 * Domain.Size;

	Domain.Origin[0] = Domain.Origin[1] = Domain.Origin[2] = 0;


	} // omp single

	return ;
}

#else // ! PERIODIC

static void set_computational_domain()
{
	find_domain_center(&Domain.Center[0]);
	
	find_largest_particle_distance(&Domain.Size);

	#pragma omp for
	for (int i = 0; i < 3; i++)
		Domain.Origin[i] = Domain.Center[i] - 0.5 * Domain.Size;

	return ;
}


/*
 * Domain Center is not the center of mass but the median of mass, which is
 * less sensitive to outliers.
 */

static Float * restrict buffer = NULL;

static void find_domain_center(double Center_Out[3])
{
	Float center[3] = { 0 };

	#pragma omp single
	buffer = Malloc(Task.Npart_Total * sizeof(*buffer), "buffer");
	
	for (int i = 0; i < 3; i++) {
		
		#pragma omp for
		for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
			buffer[ipart] = P.Pos[i][ipart];

		center[i] = Median(Task.Npart_Total, buffer);
		
		#pragma omp barrier
	}

	#pragma omp single
	buffer = Realloc(buffer, NRank * sizeof(*buffer), "buffer");

	for (int i = 0; i < 3; i++) {

		#pragma omp single
		MPI_Gather(&center[i], 1, MPI_MYFLOAT, buffer, 1, MPI_MYFLOAT,
				   Master, MPI_COMM_WORLD);

		Center_Out[i] = Median(NRank, buffer);
		
		#pragma omp barrier
	}

	#pragma omp single
	{

	MPI_Bcast(&Center_Out[0], 3, MPI_DOUBLE, Master, MPI_COMM_WORLD);
	
	Free(buffer);

	} // omp single

	return ;
}

static double Max_Distance = 0;

static void find_largest_particle_distance(double *size_out)
{
	#pragma omp single
	Max_Distance = 0;

	#pragma omp for reduction(max:Max_Distance)
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
		
		double dx = fabs(P.Pos[0][ipart] - Domain.Center[0]);
		double dy = fabs(P.Pos[1][ipart] - Domain.Center[1]);
		double dz = fabs(P.Pos[2][ipart] - Domain.Center[2]);

		Max_Distance = fmax(Max_Distance, dx);
		Max_Distance = fmax(Max_Distance, dy);
		Max_Distance = fmax(Max_Distance, dz);

	} // for ipart

	#pragma omp single
	{
	
	MPI_Allreduce(MPI_IN_PLACE, &Max_Distance, 1, MPI_DOUBLE, MPI_MAX,
		MPI_COMM_WORLD);

	*size_out = 2.001 * Max_Distance; // 2.001 helps with cancellation

	}

	return ; 
}

#endif // ! PERIODIC

#ifdef DEBUG_DOMAIN

static int compare_bunches_by_target(const void *a, const void *b)
{
	const struct Bunch_Node *x = (const struct Bunch_Node *) a;
	const struct Bunch_Node *y = (const struct Bunch_Node *) b;

	return (int) (x->Target < y->Target) - (x->Target > y->Target);
}

static void print_domain_decomposition (const int max_level)
{

	Qsort(&D[0], NBunches, sizeof(*D), &compare_bunches_by_target);

	#pragma omp master
	{

	printf(" No | Split | npart  |   sum  |   first    | trgt  | lvl |"
			"  Cost  | CumCost | Max PH key,   Max_level %d \n", max_level);

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

		printf("%3d |   %d   | %06llu | %06zu | %010d | %05d | %03d | %06g "
				"| %06g || ", i, D[i].Bunch.Modify, D[i].Bunch.Npart, sum,
				D[i].Bunch.First_Part, D[i].Bunch.Target, D[i].Bunch.Level,
				D[i].Bunch.Cost, csum);

		Print_Int_Bits64(D[i].Bunch.Key);
	}

	Assert(sum == Sim.Npart_Total, "More or less particles in D than in Sim");

	printf("\nDistribution Table : \n"
		 	"Task |   Cost   |   Imbal   |   Npart   |   Imbal   |   "
			"Split_Idx \n");

	for (int i = 0; i < NTarget; i++) 
		printf("%4d | %5.03g | %+5.03g | %9d | %+5.03g | %8d\n",
				i, Cost[i], (Cost[i]-Mean_Cost)/Mean_Cost, 
				Npart[i], ((double)Npart[i] - Mean_Npart)/Mean_Npart, 
				Split_Idx[i]);

	} // omp master

	#pragma omp barrier

	return ;
}
#endif // DEBUG_DOMAIN
