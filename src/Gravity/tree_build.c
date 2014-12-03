#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#ifdef GRAVITY_TREE

#define NODES_PER_PARTICLE 0.7 

struct Tree_Node *Tree = NULL;
int NNodes = 0;
int Max_Nodes = 0;

void gravity_tree_init();

static void build_top_tree();
static inline void create_node_from_bunch(const int, const int, 
											const peanoKey, const int, int*);

static int build_subtree(const int, const int, const int, const int);
static inline bool particle_is_inside_node(const peanoKey,const int,const int);
static inline void add_particle_to_node(const int, const int);
static inline void create_node_from_particle(const int, const int, 
											const peanoKey, const int, int*);

static void finalise_tree();
static inline void add_node_to_node(const int, const int);

static inline int key_fragment(const int node);

/*
 * This builds the tree. 
 * First the tree in every local node is build starting from the bunchnode. 
 * In the global tree memory, all bunchnodes are seperated by npart nodes.
 * Then the tree from the top nodes is constructed and the local nodes 
 * are moved in blocks.
 *
 */

void Gravity_Tree_Build()
{
	Profile("Build Gravity Tree");

	gravity_tree_init();
	// build_top_tree();
	// finalise_tree();
	
	//rprintf("Top tree has %d node, maximum depth %d \n");

	#pragma omp single
	{
	
	NNodes = build_subtree(0, Task.Npart_Total, 0, 0);

	printf("\nTree build: %zu nodes for %d particles\n", 
			NNodes, Task.Npart_Total );

	} // omp single
	
	finalise_tree();
	
	Sig.Force_Tree_Build = false;

	Profile("Build Gravity Tree");

	return ;
}

static void build_top_tree()
{	
	NNodes = 0;

	create_node_from_bunch(0, 0, 0, 0, &NNodes);

	Tree[0].Pos[0] = Domain.Center[0];
	Tree[0].Pos[1] = Domain.Center[1];
	Tree[0].Pos[2] = Domain.Center[2];

	int last_parent = 0;

	uint64_t last_key = B[0].Key >> 3;

	#pragma omp single nowait
	for (int i = 1; i < NBunches; i++) {
	
		int node = 0;
		int lvl = 0;
		int parent = node;

		uint64_t key = B[i].Key;

		bool bunch_starts_new_branch = true;

		for (;;) {
			
			if (particle_is_inside_node(key, lvl, node)) {
			
				if (Tree[node].DNext < 0) { // points to a rank -> leaf
				
					Tree[node].DNext = 0;

					create_node_from_bunch(i-1,node,last_key,lvl+1,&NNodes);

					last_key >>= 3;
				}

				add_node_to_node(i, node);
				
				bunch_starts_new_branch &= (node != last_parent);

				node++;
				lvl++;
				key >>= 3;
			
			} else { // skip
	
				if (Tree[node].DNext == 0 || node == NNodes - 1)   
					break; // reached end of branch or top tree
				
				node += fmax(1, Tree[node].DNext);
			}
		} // for (;;)

		if (bunch_starts_new_branch && (B[i].Target < 0) ) { // kick off ?
		
			#pragma omp task
			build_subtree(-(B[i].Target+1), B[i].Npart, parent, parent);

			NNodes += NODES_PER_PARTICLE * B[i].Npart; // space for subtree

			Tree[parent].DNext = NODES_PER_PARTICLE * B[i].Npart;
		}

		if (Tree[node].DNext == 0) 				// set DNext for internal node
			Tree[node].DNext = NNodes - node; 	// only delta
	
		create_node_from_bunch(i, parent, key, lvl, &NNodes); // sibling
	
		last_key = key >> 3;
		last_parent = parent;

	} // for (i < NBunches)

	#pragma omp barrier
	
	return ;
}



/* 
 * We complete the tree by walking the topnode tree backwards adding bunch node
 * properties. 
 * After the initial build, some inner DNext pointers are 0. This is corrected 
 * setting these pointers and closing the P-H curve through the tree.
 * Finally we do some final operations on the node contents. Computation can 
 * be nicely overlapped with Open MP.
 */

static void finalise_tree()
{
	//#pragma omp single nowait // close PH loop
	{

	Tree[0].DNext = 0; 

	int stack[N_PEANO_TRIPLETS + 1] = { 0 }; 
	int lowest = 0;

	for (int i = 1; i < NNodes; i++) {
		
		int lvl = Level(i);

		while (lvl <= lowest) { // set pointers

			int node = stack[lowest];
	
			if (node > 0)
				Tree[node].DNext = i - node;

			stack[lowest] = 0;

			lowest--;
		} 
		
		if (Tree[i].DNext == 0) { // add node to stack
			
			stack[lvl] = i;
			
			lowest = lvl;
		}
		
	} // for
	
	} // omp single nowait

	#pragma omp for nowait
	for (int i = 0; i < NBunches; i++) {
		
		int src = B[i].Target;

		if (src < 0)
			continue;
	
		int node = src;

		while (node != 0) {
		
			node += Tree[node].DUp;

			#pragma omp critical
			add_node_to_node(i, node);
		
		} // while
	} // for i

	#pragma omp for
	for (int i = 0; i < NNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;
	}

	return ;
}



/*
 * A subtree is build starting at node index "top". We use that the particles 
 * are in Peano-Hilbert order, i.e. that a particle will branch off as late as 
 * possible from the previous one.
 * In the tree, DNext is the difference to the next sibling in the walk, if the
 * node is not opened. Opening a node is then node++. If DNext is negative, it 
 * points to Npart particles starting at ipart=-DNext-1, and the next node in 
 * line is node++. DNext=0 is only once per level, at the end of the branch. 
 * The tree saves only one particle per node, up to eight are combined in a 
 * node. This is achieved on the fly in an explicit cleaning step when a
 * particle opens a new branch. 
 * If the tree reaches the dynamic range of the 128bit Peano-Hilbert key, it 
 * dumps all particles in the level 42 node, making the algorithm effectively 
 * N^2  again. This happens at a depth, which corresponds to a distance of 
 * Domain.Size/2^42, hence only occurs with double precision positions. The 
 * Tree.Bitfield contains the level of the node and the Peano-Triplet of the 
 * node at that level. See Tree definition in gravity.h. 
 */

static int build_subtree(const int istart, const int npart, const int offset, 
		const int top)
{
#ifdef DEBUG
	printf("DEBUG: (%d:%d) Sub-Tree Build istart=%d npart=%d offs=%d top=%d \n"
			, Task.Rank, Task.Thread_ID, istart, npart, offset, top);
#endif

	int nNodes = 0; // local in this subtree
	
	create_node_from_particle(istart, offset, 0, top, &nNodes); // first 

	Node_Set(TOP,0);
	
	Tree[0].Pos[0] = Domain.Center[0];
	Tree[0].Pos[1] = Domain.Center[1];
	Tree[0].Pos[2] = Domain.Center[2];

	int last_parent = offset;		// last parent of last particle

	Float px = (P[0].Pos[0] - Domain.Origin[0]) / Domain.Size;
	Float py = (P[0].Pos[1] - Domain.Origin[1]) / Domain.Size;
	Float pz = (P[0].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
	peanoKey last_key = Reversed_Peano_Key(px, py, pz);
	
	last_key >>= 3; 

	for (int ipart = istart+1; ipart < Task.Npart_Total; ipart++) {
		
		double px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		double py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		double pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		peanoKey key = Reversed_Peano_Key(px, py, pz);

		int node = offset; //+1		// current node
		int lvl = top;	//+1		// counts current level
		int parent = node;			// parent of current node

		bool ipart_starts_new_branch = true; // flag to remove leaf nodes
		
		while (lvl < N_PEANO_TRIPLETS) {
			
			if (particle_is_inside_node(key, lvl, node)) { // open node	
				
				if (Tree[node].Npart == 1) { 	// refine 
	
					Tree[node].DNext = 0;		

					create_node_from_particle(ipart-1, node, last_key, lvl+1, 
							&nNodes); // son of node

					last_key >>= 3;
				}  
				
				add_particle_to_node(ipart, node); // add ipart to node

				ipart_starts_new_branch &= (node != last_parent);

				parent = node;

				node++; // decline into node
				lvl++;
				key >>= 3;

			} else { // skip node
				
				if (Tree[node].DNext == 0 || node == nNodes - 1)   
					break; // reached end of branch
				
				node += fmax(1, Tree[node].DNext);
			}
		} // while (lvl < 42)

		if (lvl > N_PEANO_TRIPLETS-1) {	// particles closer than PH resolution
		
			P[ipart].Tree_Parent = parent;
			
			continue; 					// tree cannot be deeper
		}
		
		if (ipart_starts_new_branch) { // collapse particle leaf nodes

			int n = 0; 

			if (Tree[node].Npart <= 8) 
				n = node;			
			else if (Tree[last_parent].Npart <= 8) 
				n = last_parent;
			
			if (n != 0) {

				Tree[n].DNext = -ipart + Tree[n].Npart - 1;
					
				int nZero = nNodes - n - 1;

				nNodes = n + 1;

				memset(&Tree[nNodes], 0, nZero*sizeof(*Tree));

				int first = -(Tree[n].DNext + 1); // correct parent pointer
				int last = first + Tree[n].Npart;

				for (int jpart = first; jpart < last; jpart++)
					P[jpart].Tree_Parent = n;
			}	
		}
	
		if (Tree[node].DNext == 0) 				// set DNext for internal node
			Tree[node].DNext = nNodes - node; 	// only delta
			
		create_node_from_particle(ipart, parent, key, lvl, &nNodes); // sibling
	
		last_key = key >> 3;
		last_parent = parent;
	
	} // for ipart

	return nNodes;
}

/*
 * For particle and node to overlap the peano key triplet at this tree level 
 * has to be equal. Hence the tree cannot be deeper than the PH key
 * resolution, which for our 128 bits length is 42. This corresponds to 
 * distances less than 2^-42, small enough for single precision.
 */

static inline bool particle_is_inside_node(const peanoKey key, const int lvl,		const int node)
{
	int part_triplet = key & 0x7;

	int node_triplet = key_fragment(node); 

	return (node_triplet == part_triplet); 
}

/*
 * For top nodes negative DNext points to other MPI Ranks, so during the walk
 * the particle has to be exported there.
 */

static inline void create_node_from_bunch(const int i, const int parent, 
		const peanoKey key, const int lvl, int *nNodes)
{
	const int node = (*nNodes)++;

	int keyfragment = (key & 0x7) << 6;

	Tree[node].Bitfield = lvl | keyfragment;

	Node_Set(TOP, node); // bunches always make top nodes

	if (B[i].Target == Task.Rank) 
		Node_Set(LOCAL, node);
	 else 
		Tree[node].DNext = -B[i].Target; // negative DNext goes to MPI Rank

	const int sign[3] = { -1 + 2 * (B[i].Pos[0] > Tree[parent].Pos[0]),
	 			     	  -1 + 2 * (B[i].Pos[1] > Tree[parent].Pos[1]),
	 			          -1 + 2 * (B[i].Pos[2] > Tree[parent].Pos[2]) }; 
	
	Float size = Domain.Size / (1ULL << lvl);

	Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size * 0.5;
	Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size * 0.5;
	Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size * 0.5;

	Tree[node].DUp = node - parent;

	return ;
}

/*
 * We always add nodes at the end of the subtree. Particle pointers are 
 * negative and offset by one to leave DNext=0 indicating unset.
 */

static inline void create_node_from_particle(const int ipart,const int parent, 
		const peanoKey key, const int lvl, int *nNodes)
{
	const int node = (*nNodes)++;

#ifdef DEBUG
	Assert(node < Max_Nodes, "Too many nodes (%d>%d), increase "
			"NODES_PER_PARTICLE=%g", nNodes, Max_Nodes, NODES_PER_PARTICLE);
#endif

	Tree[node].DNext = -ipart - 1;

	int keyfragment = (key & 0x7) << 6;

	Tree[node].Bitfield = lvl | keyfragment | (1UL << 9);

	const int sign[3] = { -1 + 2 * (P[ipart].Pos[0] > Tree[parent].Pos[0]),
	 			     	  -1 + 2 * (P[ipart].Pos[1] > Tree[parent].Pos[1]),
	 			          -1 + 2 * (P[ipart].Pos[2] > Tree[parent].Pos[2]) }; 
	
	Float size = Domain.Size / (1 << lvl);

	Tree[node].Pos[0] = Tree[parent].Pos[0] + sign[0] * size * 0.5;
	Tree[node].Pos[1] = Tree[parent].Pos[1] + sign[1] * size * 0.5;
	Tree[node].Pos[2] = Tree[parent].Pos[2] + sign[2] * size * 0.5;

	Tree[node].DUp = node - parent;

	P[ipart].Tree_Parent = node;

	add_particle_to_node(ipart, node); 

	return ;
}

/*
 * Adding a particle or a node to a node is really the same thing ...
 */

static inline void add_particle_to_node(const int ipart, const int node)
{
	Tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	Tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	Tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;

	Tree[node].Mass += P[ipart].Mass;

	Tree[node].Npart++;
	
	return ;
}

static inline void add_node_to_node(const int src, const int target)
{
	Tree[target].CoM[0] += Tree[src].CoM[0];
	Tree[target].CoM[1] += Tree[src].CoM[1];
	Tree[target].CoM[2] += Tree[src].CoM[2];

	Tree[target].Mass += Tree[src].Mass;

	Tree[target].Npart += Tree[src].Npart;

	return ;
}


/*
 * These functions set/extract bits from the bitmask to control
 * the tree construction and walk. They will be inlined by the compiler and
 * cost only a few cycles.  
 */

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (Tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

int Level(const int node)
{
	return Tree[node].Bitfield & 0x3FUL; // return but 0-5
}

bool Node_Is(const enum Tree_Bitfield bit, const int node)
{
	return Tree[node].Bitfield & (1UL << bit);
}

void Node_Set(const enum Tree_Bitfield bit, const int node)
{
	Tree[node].Bitfield |= 1UL << bit;

	return ;
}

void Node_Clear(const enum Tree_Bitfield bit, const int node)
{
	Tree[node].Bitfield &= ~(1UL << bit);

	return ;
}


/*
 * Initialise the first tree node
 */

void gravity_tree_init()
{
	#pragma omp single
	{

	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	size_t nBytes = Max_Nodes * sizeof(*Tree);

	if (Tree == NULL)
		Tree = Malloc(nBytes, "Tree");
	
	memset(Tree, 0, nBytes);

	} // omp single

	return ;
}



void test_gravity_tree(const int nNodes)
{
	for (int node = 0; node < nNodes; node++) {
	
		int lvl = Level(node);
	
		int n = node + 1;

		float mass = 0;
		int npart = 0;
		int nout = 0;

		double nSize = Domain.Size / (float)(1ULL << lvl);

		if (Tree[node].DNext < 0)
			continue;

		while (Level(n) > lvl) {

			if (Tree[n].DNext < 0) {
			
				int first = -Tree[n].DNext - 1;
				int last = first + Tree[n].Npart;

				for (int jpart = first; jpart < last; jpart++ ) {

				 	npart++;

					mass += P[jpart].Mass;
	
					float dx = fabs(P[jpart].Pos[0] - Tree[n].Pos[0]);
					float dy = fabs(P[jpart].Pos[1] - Tree[n].Pos[1]);
					float dz = fabs(P[jpart].Pos[2] - Tree[n].Pos[2]);

					if (dx > nSize * 0.5) 
						if (dy > nSize * 0.5) 
							if (dz > nSize * 0.5)
								nout++;
				}
			}

			n++; 
		}

		printf("%d m=%g,%g N=%d,%d nsize=%g nout=%d\n", node, mass, 
				Tree[node].Mass, npart, Tree[node].Npart, 	nSize, nout);
	}

	return ;
}


#endif // GRAVITY_TREE
