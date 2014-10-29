#include "../globals.h"
#include "../domain.h"
#include "gravity.h"

#include "memory.h"
#define NODES_PER_PARTICLE 1.1 // Springel 2005

struct Tree_Node {
	uint32_t Bitfield; 
	// bit 0-5:level, 6-8:key, 9-15:nLeaves, 16: DNext dangling, 17-32:free
	int DNext;		 // this is the DELTA to the next node; or -ipart
	Float CoM[3];
	float Mass;
	int Npart;
	int DUp;
} *Tree = NULL;

static size_t NNodes = 0;
static size_t Max_Nodes = 0;

static inline int level(const int node); // bitmask functions
static inline int key_fragment(const int node);
static inline Float node_size(const int node);
static inline bool dNext_is_dangling(const int node);
static inline void set_dNext_dangling(const int node);
static inline void unset_dNext_dangling(const int node);
static inline void add_particle_to_bitfield(const int node);
static inline int nleafs_in_node(const int node); 

static inline void add_node(const int ipart, const int node);
static inline bool particle_is_inside_node(const int ipart, const int node);
static inline void add_particle_to_node(const int ipart, const int node);

static void print_int_bits64(const uint64_t val)
{
	for (int i = 63; i >= 0; i--) {
		printf("%llu", (val & (1ULL << i) ) >> i);
		if (i % 3 == 0 && i != 0)
			printf(".");
	}
	printf("\n");fflush(stdout);

	return ;
}

/*
 * This builds the tree. 
 */

void Build_Tree()
{
	Profile("Build Gravity Tree");

	Init_Tree();

	/* Build_Top_Tree();
	
	   for (int i = 0; i < NBunches; i++) {

	   int  node_start, ipart_start, ipart_stop;
		
	*/
 
 	/* Here the top tree is complete, now build local nodes downwards */

for (int n=0; n<NNodes; n++)
	printf("START n=%d  np=%d next=%d up=%d mass=%g \n", 
			n,  Tree[n].Npart, Tree[n].DNext, Tree[n].DUp,Tree[n].Mass);
bool tmp = particle_is_inside_node(0, 0);


	add_node(0, 0); // add the first particle by hand
	Tree[0].Bitfield &= ~0x1FF;

	for (int ipart = 1; ipart < 20; ipart++) {

printf("\n\nIPART=%d NNodes=%zu \n", ipart, NNodes);

		int node = 0;  			// current node
		int lvl = 0;			// counts current level
		int parent = 0;		// parent of current node
		int key = 0;			// peano key of ipart
int cnt = 0;
		for (;;) {

			int part_triplet = (key >> lvl) & 0x7; 
//if (cnt ++ > 10) break;

//printf("\nLOOP : node=%d npart=%d next=%d level=%d dang=%d\n", node, Tree[node].Npart, Tree[node].DNext, level(node), dNext_is_dangling(node));

if (level(node) != lvl)
printf("ERROR LEVEL! %d %d \n", level(node), lvl);

			if (particle_is_inside_node(ipart, node)) { // open node
			
//printf("IN : \n");
	
				if (Tree[node].Npart == 1) { // refine
				
					int jpart = -Tree[node].DNext; // save current particle
	
//printf("REFINE : node=%d npart=%d jpart=%d dangle=%d\n", 
//node,Tree[node].Npart, jpart, dNext_is_dangling(node));
					
					set_dNext_dangling(node); // DNext is now free
					
					add_node(jpart, node); // add  daughter to current node
				}  
				
				add_particle_to_node(ipart, node); // ipart to current node

//printf("DECLINE %d->%d  \n", node, node+1);
				parent = node;
				node++; // decline in node containing jpart
				lvl++;

			} else { // skip

//printf("OUT %d : \n", Tree[node].DNext);

				if (dNext_is_dangling(node) || node == NNodes - 1) {  

//printf("LEAF : %d \n", dNext_is_dangling(node));

					break; // reached end of my branch 
				} 

				if (Tree[node].Npart == 1) { // internal leaf

					node++; // might need to split its sibling 
//printf("NEXT : +%1\n");
				} else {
			
//printf("NEXT : +%d\n",Tree[node].DNext);
					node += Tree[node].DNext; // skip whole subtree

				}
			}

		} // for (;;)

	//	if (parent != last_parent) { // compactify tree
//printf("COMPACTIFY node=%d, dNext=%d NNodes=%d newdNext=%d \n",
//node, Tree[node].DNext, NNodes, NNodes-node);
	//		NNodes = parent; // all trees smaller than 8 are particles
	//		}
	
		if (dNext_is_dangling(node)) { // set internal dangling next node

			Tree[node].DNext = NNodes - node; // only delta

			unset_dNext_dangling(node);
		}
		//last_parent = parent;
		int parent2 = node - Tree[node].DUp;

	printf("%d %d \n", parent, parent2);
//printf("ADD OUT : node=%d npart=%d parent=%d lvl=%d, dangling dNext=%d\n", 
//node, Tree[node].Npart, parent, level(node), dNext_is_dangling(node));

		add_node(ipart, parent); // add a sibling to current node

	} // for ipart


for (int n=0; n<NNodes; n++) {
	//if (Tree[n].DNext != -1 && dNext_is_dangling(n)) {
printf("TEST n=%d np=%d next=%d up=%d mass=%g level=%d dangle=%d \n", 
n,  Tree[n].Npart, Tree[n].DNext, Tree[n].DUp,Tree[n].Mass, level(n), 
dNext_is_dangling(n));
print_int_bits64(Tree[n].Bitfield);
//	}
}
printf("\n");
for (int ipart = 0; ipart < 20; ipart++) 
print_int_bits64(P[ipart].Peanokey);

	Profile("Build Gravity Tree");
Print_Memory_Usage();
exit(0);
	return ;
}

/*
 * For particle and node to overlap the peano key triplet at this
 * tree level has to be equal. If the tree is deeper than the number of 
 * bit triplets in the peano key (i.e. level>21), we always return 
 * "false", causing the particles to be added as leaves at the end of the 
 * tree in random order. The tree walk then has to look at all of them
 * which make this N^2.
 */

static inline bool particle_is_inside_node(const int ipart, const int node)
{
	const int lvl = level(node);

	const int shift = 63 - 3 * lvl; // bits 0-5

	const uint64_t part_mask = 0x07ULL << shift; // shift first 3 bits

	uint64_t part_triplet = (P[ipart].Peanokey & part_mask) >> shift;

	uint64_t node_triplet = key_fragment(node); 

	return (node_triplet == part_triplet) && (lvl < 22); // branch free
}

/*
 * We always add nodes at the end of the tree
 */

static inline void add_node(const int ipart, const int parent)
{
	const int node = NNodes++;

	unset_dNext_dangling(node);

	Tree[node].DNext = -ipart;
	
	uint32_t lvl = level(parent) + 1; // construct bitfield
	int shift = 63 - 3*lvl;

	uint64_t tmp_key = (P[ipart].Peanokey & (0x07ull << shift));
	uint32_t keyfragment = (uint32_t) (tmp_key >> (shift - 6)); 
	
	Tree[node].Bitfield = lvl | keyfragment; 

	Tree[node].DUp = node - parent;

	add_particle_to_node(ipart, node); 

//printf("ADD ipart=%d n=%d parent=%d lvl=%d shft=%d Par_UP=%d Pdangle=%d \n", 
//ipart, node, parent,lvl, shift, Tree[parent].DUp, dNext_is_dangling(parent));
//print_int_bits64(Tree[node].Bitfield);
	return ;
}

static inline void add_particle_to_node(const int ipart, const int node)
{
	add_particle_to_bitfield(node); 

	Tree[node].CoM[0] += P[ipart].Pos[0] * P[ipart].Mass;
	Tree[node].CoM[1] += P[ipart].Pos[1] * P[ipart].Mass;
	Tree[node].CoM[2] += P[ipart].Pos[2] * P[ipart].Mass;
	
	Tree[node].Mass += P[ipart].Mass;

	/* abuse the node below */

	Tree[node].Npart++;

	return ;
}

/*
 * These functions set/extract bits from the bitmask to control
 * the tree construction and walk. They will be inlined by the compiler and
 * cost only a few cycles. We store bitmasks separately to ensure the
 * correct type of the constant.
 */

static inline int level(const int node)
{
	const uint32_t bitmask = 0x3F;

	return Tree[node].Bitfield & bitmask; // return but 0-5
}

static inline int key_fragment(const int node)
{
	const uint32_t bitmask = 7UL << 6;

	return (Tree[node].Bitfield & bitmask) >> 6; // return bit 6-8
}

static inline int nleafs_in_node(const int node) // count up to 2^6-1
{
	const uint32_t bitmask = 15UL << 9;

	return (Tree[node].Bitfield & bitmask) >> 9; // return bit 9-15
}

static inline void add_particle_to_bitfield(const int node)
{
	uint32_t nleafs = 1 + nleafs_in_node(node);

	const uint32_t bitmask = 15UL << 9;
	Tree[node].Bitfield &= ~bitmask; // clear bits 9-15

	Tree[node].Bitfield |= nleafs << 9; // set new bits 9-15

	return ;
}

static inline bool dNext_is_dangling(const int node)
{
	return (Tree[node].Bitfield & (1UL << 15)) >> 15; // return bit 16
}

static inline void set_dNext_dangling(const int node)
{
	//Tree[node].DNext = 0xFFFFFFF;
	Tree[node].Bitfield |= 1UL << 15; // set bit 16 = 1
	
	return ; 
}

static inline void unset_dNext_dangling(const int node)
{
	Tree[node].Bitfield &= ~(1UL << 15); // clear bit 16 
	
	return ; 
}


static inline Float node_size(const int node)
{
	int lvl = level(node);

	return Domain.Size / ((Float) (1UL << lvl)); // Domain.Size/2^level
}


void Init_Tree()
{
	#pragma omp single
	{

	Max_Nodes = Task.Npart_Total_Max * NODES_PER_PARTICLE;
	
	size_t nBytes = Max_Nodes * sizeof(*Tree);

	if (Tree == NULL)
		Tree = Malloc(nBytes, "Tree");
	
	memset(Tree, 0, nBytes);

	Tree[0].CoM[0] = 0; // will hold global CoM of sim
	Tree[0].CoM[1] = 0;
	Tree[0].CoM[2] = 0;

	Tree[0].Npart = 0;
	Tree[0].DUp = -1;
	Tree[0].Mass = 0;
	Tree[0].DNext = -1;
	Tree[0].Bitfield = 0;

	NNodes = 0;

	} // omp single

	return ;
}

/*
 * This function computes the gravitational acceleration by walking the tree  
 */

void Accel_Gravity_Tree(const int ipart, double *acc, double *pot)
{

	if (Sig.Fullstep) // need to build new one
		Free(Tree);

	return ;
}
