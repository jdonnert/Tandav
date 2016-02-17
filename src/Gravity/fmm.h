/*
 * Gravity using the Fast Multipole Method with dual tree traversal.
 */

#if defined(GRAVITY) && defined(GRAVITY_TREE)

extern struct FMM_Nodes {
	int * restrict DNext;
	uint32_t * restrict Bitfield;
	int * restrict DUp;
	int * restrict Npart;
	Float * restrict Pos[3];
	Float * restrict Mass;
	Float * restrict CoM[3];
	Float * restrict Dp[3];
} Nodes;

uint32_t NNodes;

extern struct FMM_Leafs { // vectors of particles and a parent node
	int * restrict First;
	int * restrict N;
	int * restrict Up;
} Leafs;

uint32_t NLeafs;

#endif // GRAVITY && GRAVITY_FMM
