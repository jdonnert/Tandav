#include "globals.h"
#include "domain.h"
#include "peano.h"

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

 	Sort_Particles_By_Peano_Key();

	for (;;) {
	
		fill_bunchlist();

		double imbalance = compute_imbalance();

		if ( imbalance < 0.05 )
			break;
		
		split_bunches();
	} 
	
	communicate_particles();

	Profile("Domain Decomposition");

	return ;
}

void Init_Domain_Decomposition()
{
	int nTask = Sim.NTask, nThreads = Sim.NThreads;
	
	size_t nBytes = (1ULL << DOMAIN_CEILING) * Sim.NTask * sizeof(Bunchlist); 

	Bunchlist = Malloc(nBytes);

	memset(Bunchlist, 0, nBytes);

	/* Init the first bunches */

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
