#include "globals.h"
#include "proto.h"

static const int nTask = Sim.NTask, nThreads = Sim.NThreads;

typedef bunch struct bunch_information {
	int npart[NPARTYPES];
	int next_bunch; 
	int target_task;
	peanokey key_max;
	peanokey key_min;
};

/* This function distributes particles in bunches, which are continuous
 * on the Peano curve. Because of this the bunches correspond to parts
 * of the tree and adjacent bunches have ghost nodes in the tree. 
 * We keep a global list of all the bunches that also contains the workload and 
 * memory footprint of each bunch. An optimal way of distributing bunches
 * minimises memory and workload imbalance over all Tasks. To approach this
 * we keep more bunches than Tasks */
void Domain_Decomposition()
{
	bunch *bunchlist = NULL;

	// construct bunches
	// 		one bunch per task
	//
	// 		do  
	//			construct mean mem footprint & workload
	//			
	//			split bunches over threshold in 8 (go down the tree on step)
	//		while (n_above_threshold)	
	//
	// sort particles on tasks
	// export ghostlist for tree construction

	return ;
}
