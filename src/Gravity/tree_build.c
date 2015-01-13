#include "../globals.h"
#include "gravity.h"
#include "../domain.h"

#ifdef GRAVITY_TREE

#define NODES_PER_PARTICLE 0.5 

struct Tree_Node *Tree = NULL;
int NNodes = 0;
int NTop_Nodes = 0;
int Max_Nodes = 0;

void gravity_tree_init();

static int transform_bunch_into_top_node(const int);

static void build_subtree(const int, const int, const int, int *);
static void finalise_subtree(const int, const int, const int);

static inline bool particle_is_inside_node(const peanoKey,const int,const int);
static inline void add_particle_to_node(const int, const int);
static inline void create_node_from_particle(const int, const int, 
		const peanoKey, const int, const int);
static peanoKey create_first_subtree_node(const int, const int, const int);

static inline int key_fragment(const int);

/*
 * This builds the tree in parallel.
 * First every local bunch is converted into a topnode and broadcasted. 
 * If the number of particles in that top node is <= 8  the topnode target 
 * points directly to the particles and the tree walk will use particles 
 * directly from the topnode. Else a subtree is build starting from the top 
 * node. Every subtree allocated a fixed amount of memory. I.e. in the global 
 * tree memory, all bunch sub trees are separated by 
 * TNode.Npart * NODES_PER_PARTICLE nodes.
 */

void Gravity_Tree_Build()
{
	Profile("Build Gravity Tree");

	gravity_tree_init();

	NTop_Nodes = NBunches;

	#pragma omp for schedule(static,1)
	for (int i = 0; i < NBunches; i++) {
	
		int level = 0, src = 0;
		printf("TopNode %d, Target=%d \n", i, D[i].TNode.Target);
		if (D[i].TNode.Target >= 0) { // local node
		
			level = transform_bunch_into_top_node(i);
			
			src = Task.Rank;

		} else { // non-local node
			
			src = -D[i].TNode.Target;
		}

	/*	MPI_Request *request = NULL;
		float *target = &D[i].TNode.Pos[0];
		int nBytes = (&D[i].TNode.Dp[2] - target) + sizeof(float);
	MPI_Ibcast(target, nBytes, MPI_BYTE, src, MPI_COMM_WORLD, request); */
	
		if (D[i].TNode.Target < 0)  // not local
			continue; 

		if (D[i].TNode.Npart <= 8) // No subtree needed, target stays ipart
			continue;
	
		int nNodes_subtree = ceil(D[i].TNode.Npart * NODES_PER_PARTICLE);

		#pragma omp critical
		{
		
		D[i].TNode.Target = NNodes;

		NNodes += nNodes_subtree;
		
		Assert(NNodes < Max_Nodes, "Too many nodes (%d>%d), increase "
			"NODES_PER_PARTICLE=%g", NNodes, Max_Nodes, NODES_PER_PARTICLE);
		
		} // omp critical


		int ipart = D[i].TNode.Target;

		int nlocal = 0;

printf("Tree Build: ipart=%d bunch=%d lvl=%d NNodes=%d subnodes=%d\n", ipart, i, level, NNodes, nNodes_subtree );

		#pragma omp task
		build_subtree(ipart, i, level, &nlocal);

		printf("Top Node %d has %d/%d tree nodes, nnodes per part = %g \n",
				i, nlocal, nNodes_subtree, (double)nlocal / D[i].TNode.Npart);
	}

	#pragma omp barrier

	Sig.Force_Tree_Build = false;

	Profile("Build Gravity Tree");

	return ;
}

/*
 * Set the "TNode" part of the D unions. From here onwards the members 
 * are to be understood as a Top Node, not a bunch.
 */

static int transform_bunch_into_top_node(const int i)
{	
	int ipart = D[i].Bunch.First_Part;

	D[i].TNode.Target = ipart; // save the particle index

	double px = P[ipart].Pos[0] - Domain.Origin[0]; // construct from particle
	double py = P[ipart].Pos[1] - Domain.Origin[1];
	double pz = P[ipart].Pos[2] - Domain.Origin[2];

	int level = D[i].Bunch.Level;

	double size = Domain.Size / (1ULL << level);

	D[i].TNode.Pos[0] = (floor(px/size) + 0.5) * size + Domain.Origin[0];
	D[i].TNode.Pos[1] = (floor(py/size) + 0.5) * size + Domain.Origin[1];
	D[i].TNode.Pos[2] = (floor(pz/size) + 0.5) * size + Domain.Origin[2];

printf("Transform: ipart=%d pr=%g %g %g lvl=%d size=%g pn=%g %g %g \n", ipart, px, py, pz, level, size, D[i].TNode.Pos[0],D[i].TNode.Pos[1],D[i].TNode.Pos[2] );

	return level;
}



/*
 * A subtree is build starting at node index "offset", from top node at index 
 * "tnode_idx". The first node is build by hand as a copy of the top node.
 * We use that the particles are in Peano-Hilbert order, i.e. that a 
 * particle will branch off as late as possible from the previous one. This 
 * means refining a node can be done via a split and reassignment of ipart and
 * ipart-1.
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

static void build_subtree(const int istart, const int tnode_idx, 
		const int top_level, int *nNodes_out)
{
//#ifdef DEBUG
	printf("DEBUG: (%d:%d) Sub-Tree Build istart=%d npart=%d offs=%d tnode=%d\n"
		,Task.Rank, Task.Thread_ID, istart, D[tnode_idx].TNode.Npart, 
		D[tnode_idx].TNode.Target, tnode_idx);
//#endif

	const int subtree_start = D[tnode_idx].TNode.Target; 

	peanoKey last_key = create_first_subtree_node(istart, tnode_idx, top_level);

	int nNodes = 1; // for this subtree

	int last_parent = subtree_start; // last parent of last particle

	last_key >>= 3; 

	for (int ipart = istart + 1; ipart < Task.Npart_Total; ipart++) {
		
		double px = (P[ipart].Pos[0] - Domain.Origin[0]) / Domain.Size;
		double py = (P[ipart].Pos[1] - Domain.Origin[1]) / Domain.Size;
		double pz = (P[ipart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
		peanoKey key = Reversed_Peano_Key(px, py, pz);

		int node = subtree_start;   // current node
		int lvl = top_level;		// counts current level
		int parent = node;			// parent of current node

		bool ipart_starts_new_branch = true; // flag to remove leaf nodes
		
		while (lvl < N_PEANO_TRIPLETS) { 
			
			if (particle_is_inside_node(key, lvl, node)) { // open node	
				
				if (Tree[node].Npart == 1) { 	// refine 
	
					Tree[node].DNext = 0;		

					int new_node = subtree_start + nNodes++;

					create_node_from_particle(ipart-1, node, last_key, lvl+1, 
							new_node); 	// new_node is a son of node

					last_key >>= 3;
				}  
				
				add_particle_to_node(ipart, node); // add ipart to node

				ipart_starts_new_branch &= (node != last_parent);

				parent = node;

				node++; // decline into node

				lvl++;
				
				key >>= 3;

			} else { // skip to next node
				
				if (Tree[node].DNext == 0 || node == nNodes - 1)   
					break; // reached end of branch
				
				node += fmax(1, Tree[node].DNext);
			}
		} // while (lvl < 42)

		if (lvl > N_PEANO_TRIPLETS-1) {	// particles closer than PH resolution
		
			P[ipart].Tree_Parent = parent;
			
			continue;
		}
		
		if (ipart_starts_new_branch) { // collapse last particle's leaf nodes

			int n = -1; 
		
			if (Tree[node].Npart <= 8)
				n = last_parent;
			else if (Tree[last_parent].Npart <= 8)
				n = node;

			if (n != -1) {

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
			
		int new_node = subtree_start + nNodes++;

		create_node_from_particle(ipart, parent, key, lvl, new_node); // sibling
	
		last_key = key >> 3;
		last_parent = parent;
	
	} // for ipart

	finalise_subtree(istart, nNodes, top_level);
	
	*nNodes_out = nNodes;

	return ;
}

/*
 * The first node in the subtree needs special treatment and will later
 * be overwritten. Its properties are identical with the topnode from the
 * domain decomposition, which makes a great consistency check.
 */

static peanoKey create_first_subtree_node(const int istart, 
		const int tnode_idx, const int top_level)
{
	const int node = D[tnode_idx].TNode.Target; // where the subtree starts

	Float px = (P[istart].Pos[0] - Domain.Origin[0]) / Domain.Size; 
	Float py = (P[istart].Pos[1] - Domain.Origin[1]) / Domain.Size; 
	Float pz = (P[istart].Pos[2] - Domain.Origin[2]) / Domain.Size;
		
	peanoKey last_key = Reversed_Peano_Key(px, py, pz) >> (3*top_level);

	create_node_from_particle(istart, node, last_key, top_level, node); 

	Tree[node].Pos[0] = D[tnode_idx].TNode.Pos[0]; // correct position
	Tree[node].Pos[1] = D[tnode_idx].TNode.Pos[1]; // because parent node 
	Tree[node].Pos[2] = D[tnode_idx].TNode.Pos[2]; // did not exist

	Tree[node].DUp = tnode_idx; // correct up pointer to lead to topnode
	
	Node_Set(TOP, node);

	return last_key;
}

/* 
 * After the build, some inner DNext pointers are 0. This is corrected 
 * setting these pointers and closing the P-H curve through the sub tree.
 * Finally we do some final operations on the node contents and remove the
 * top node of the subtree, as it is identical with D.TNode.
 */

static void finalise_subtree(const int istart, const int nNodes, 
		const int top_level)
{
	Tree[istart].DNext = 0; 

	int stack[N_PEANO_TRIPLETS + 1] = { 0 }; 
	int lowest = top_level;

	for (int i = istart + 1; i < nNodes; i++) {
		
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
	
	for (int i = istart; i < nNodes; i++) {
	
		Tree[i].CoM[0] /= Tree[i].Mass;
		Tree[i].CoM[1] /= Tree[i].Mass;
		Tree[i].CoM[2] /= Tree[i].Mass;
	}

	//copy back Mass, CoM, Dp to D[i].TNode

#ifdef DEBUG
	compare_topnode_to_subtree_head();
#endif	

	// void *src =  &Tree[subtree_start+1]; // remove topnode copy
	// void *dest = &Tree[subtree_start];
	// size_t nBytes = sizeof(*Tree) * (nNodes-1);
	// memmove(dest, src, nBytes);
	
	return ;
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
 * We always add nodes at the end of the subtree. Particle pointers are 
 * negative and offset by one to leave DNext=0 indicating unset. We assume
 * the peano key has reversed triplet order and the least significant bit 
 * carries the triplet at level "lvl".
 */

static inline void create_node_from_particle(const int ipart,const int parent, 
		const peanoKey key, const int lvl, const int node)
{
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

static inline void add_particle_to_node(const int ipart, const int node)
{
	Tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	Tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	Tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;

	Tree[node].Mass += P[ipart].Mass;

	Tree[node].Npart++;
	
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
 * Initialise tree memory
 */

void gravity_tree_init()
{
	#pragma omp single
	{

	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	NNodes = 0;
	
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
