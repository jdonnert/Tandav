#include "fmm.h"

#ifdef GRAVITY_FMM

static void prepare_fmm();
static void realloc_nodes(const int N);

uint64_t NNodes = 0, Max_Nodes = 0;
extern struct FMM_Node FMM = { NULL };
static struct FMM_Node fmm = { NULL };
#pragma omp threadprivate(fmm)
static omp_lock_t Node_Lock;
static Leaf_Nodes[] = NULL;

/*
 * This builds the FMM tree in parallel. We OpenMP decompose along top-nodes
 * and vectorize the leafs.
 */

void Gravity_FMM_Build()
{
	Profile("Grav FMM Build");

	#pragma omp single
	NNodes = 0;

	#pragma omp single
	prepare_fmm();

	fmm = FMM;

	int first_vec = D[0].TNode.First_Vec;
	int top_level = D[0].TNode.Level;

	peanoKey last_key = create_first_node(Vec[first_vec], 0, D[0].TNode.Level);
		
	NNodes = 1;

	int last_parent = 0;

	for (int i = first_vec+1; i < D[i].TNode.NVec; i++) {
	
		int ipart = Vec[i];

		peanoKey key = P.Key[ipart];

		key >>= 3 * top_level;

		int node = 0;        // running node
		int lvl = top_level; // counts current level
		int parent = node;   // parent of current node

		while (lvl < N_PEANO_TRIPLETS) {

			if (particle_is_inside_node(key, lvl, node)) { // open node	

				if (fmm[node].Npart == 1) { // refine
				
					int jpart = -fmm[node].DNext -1;

					fmm[node].DNext = 0;

					int new_node = nNodes; // is a son of "node"

					create_node_from_particle(jpart, node, last_key, lvl+1,
																new_node);
					nNodes++;

					last_key >>= 3;
				}

				parent = node;

				node++; // decline
				lvl++;
				key >>= 3;

			} else { // skip right to next node
	
				if (fmm[node].DNext == 0 || node == nNodes - 1)
					break; // reached end of branch

				node += fmax(1, tree[node].DNext);
			}
		} // while lvl

		if (fmm[node].DNext == 0)	// set DNext for internal node
			fmm[node].DNext = nNodes - node;	// only delta

		int new_node = nNodes; // is a sibling of "node"

		create_node_from_particle(ipart, parent, key, lvl, new_node);

		nNodes++;

		last_key = key >> 3;

		last_parent = parent;
			
	} // for i
	
	Profile("Grav FMM Build");

	return ;
}


void Gravity_FMM_Setup()
{
	omp_init_lock(Node_Lock);

	Max_Nodes = 0.3 * Task.Npart_Total;

	realloc_nodes(Max_Nodes, FMM);

	Leaf_Nodes = Malloc(Task.Npart_Total * sizeof(*Leaf_Node), "Leaf Nodes");

	for (int i = 0; i < NPARTYPE; i++) { // Plummer eqiv. softening
	
		Epsilon[i] = -41.0/32.0 * Param.Grav_Softening[i]; // for Dehnen K1
		Epsilon2[i] = Epsilon[i] * Epsilon[i];
		Epsilon3[i] = Epsilon[i] * Epsilon[i] * Epsilon[i];
	}

	return ;
}

void Free_Gravity_FMM(struct FMM_Node *f)
{
	Free(f->DNext);
	Free(f->Bitfield);
	Free(f->DUp);
	Free(f->Npart);
	Free(f->Pos[0]);
	Free(f->Pos[1]);
	Free(f->Pos[2]);
	Free(f->Mass);
	Free(f->CoM[0]);
	Free(f->CoM[1]);
	Free(f->CoM[2]);
	Free(f->Dp[0]);
	Free(f->Dp[1]);
	Free(f->Dp[2]);

	return ;
}

static void realloc_nodes(const int N, struct FMM_Node *f)
{
	if (f.DNext != NULL) 
		Free_Gravity_FMM(f);

	f->DNext = Malloc(N*sizeof(*FMM.DNext), "FMM.DNext");
	f->Bitfield = Malloc(N*sizeof(*FMM.Bitfield), "FMM.Bitfield");
	f->DUp = Malloc(N*sizeof(*FMM.DUp), "FMM.DUp");
	f->Npart = Malloc(N*sizeof(*FMM.Npart), "FMM.Npart");
	f->Pos[0] = Malloc(N*sizeof(*FMM.Pos[0]), "FMM.Pos[0]");
	f->Pos[1] = Malloc(N*sizeof(*FMM.Pos[1]), "FMM.Pos[1]");
	f->Pos[2] = Malloc(N*sizeof(*FMM.Pos[2]), "FMM.Pos[2]");
	f->Mass = Malloc(N*sizeof(*FMM.Mass), "FMM.Mass");
	f->CoM[0] = Malloc(N*sizeof(*FMM.CoM[0]), "FMM.CoM[0]");
	f->CoM[1] = Malloc(N*sizeof(*FMM.CoM[1]), "FMM.CoM[1]");
	f->CoM[2] = Malloc(N*sizeof(*FMM.CoM[2]), "FMM.CoM[2]");
	f->Dp[0] = Malloc(N*sizeof(*FMM.Dp[0]), "FMM.Dp[0]");
	f->Dp[1] = Malloc(N*sizeof(*FMM.Dp[1]), "FMM.Dp[1]");
	f->Dp[2] = Malloc(N*sizeof(*FMM.Dp[2]), "FMM.Dp[2]");
	
	return ;
}

static void prepare_fmm()
{
	if (NNodes != 0) { // tree build aborted, increase mem

		Max_Nodes = TREE_ENLARGEMENT_FACTOR * Max_Nodes;
		
		printf("(%d:%d) Increased tree memory to %6.1f MB, "
			"max %10d nodes, ratio %4g \n", Task.Rank, Task.Thread_ID, 
			Max_Nodes * sizeof(*Tree)/1024.0/1024.0, Max_Nodes, 
			(double) Max_Nodes/Task.Npart_Total); 
	}

	realloc_nodes(Max_Nodes, FMM);
		
	NNodes = 0;

	return ;
}

/*
 * The first node in the subtree needs special treatment. Its properties are 
 * identical with the topnode from the domain decomposition.
 */

static peanoKey create_first_node(const int first_part,
		const int tnode_idx, const int top_level)
{
	peanoKey key = P.Key[first_part];

	key >>= 3 * top_level;

	create_node_from_particle(first_part, 0, key, top_level, 0);

	fmm[0].Pos[0] = D[tnode_idx].TNode.Pos[0]; // get top node position 
	fmm[0].Pos[1] = D[tnode_idx].TNode.Pos[1]; // because parent node 
	fmm[0].Pos[2] = D[tnode_idx].TNode.Pos[2]; // did not exist

	fmm[0].DUp = tnode_idx; // correct up pointer to lead to topnode

	node_set(TOP, 0);

	return key >> 3;
}


/*
 * For particle and node to overlap the peano key triplet at this tree level 
 * has to be equal. Hence the tree cannot be deeper than the PH key
 * resolution, which for our 128 bits length is 42. This corresponds to 
 * distances less than 2^-42, small enough for single precision.
 */

static inline bool particle_is_inside_node(const peanoKey key, const int lvl,
										   const int node)
{
	int part_triplet = key & 0x7;

	int node_triplet = key_fragment(node);

	return node_triplet == part_triplet;
}

/*
 * We always add nodes at the end of the subtree. Particles are named 
 * negative and offset by one to leave DNext=0 indicating unset. We assume
 * the peano key has reversed triplet order and the least significant 3 bits 
 * carry the triplet at level "lvl".
 */

static inline void create_node_from_particle(const int ipart,const int parent,
											 const peanoKey key, const int lvl,
											 const int node)
{
	fmm[node].DNext = -ipart - 1;

	int keyfragment = (key & 0x7) << 6;

	fmm[node].Bitfield = lvl | keyfragment | (1UL << 9);

	const int sign[3] = { -1 + 2 * (P.Pos[0][ipart] > tree[parent].Pos[0]),
					      -1 + 2 * (P.Pos[1][ipart] > tree[parent].Pos[1]),
		                  -1 + 2 * (P.Pos[2][ipart] > tree[parent].Pos[2]) };

	Float size = Domain.Size / (1 << lvl);

	fmm[node].Pos[0] = fmm[parent].Pos[0] + sign[0] * size * 0.5;
	fmm[node].Pos[1] = fmm[parent].Pos[1] + sign[1] * size * 0.5;
	fmm[node].Pos[2] = fmm[parent].Pos[2] + sign[2] * size * 0.5;

	fmm[node].DUp = node - parent;

	P.Tree_Parent[ipart] = node;

	return ;
}

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

static inline void node_set(const enum Tree_Bitfield bit, const int node)
{
	tree[node].Bitfield |= 1UL << bit;

	return ;
}



#endif // GRAVITY_FMM
