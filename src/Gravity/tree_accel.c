#include "../globals.h"
#include "gravity.h"
#include "../domain.h"

#ifdef GRAVITY_TREE

#define NODE_OPEN_PARAM 0.01 // open criterion parameter

static struct GravityDataForExport {
	int Target_Task;
	Float Pos[3];
	Float Acc;
} DataOut;

static struct GravityDataForImport {
	int Target_Part;
	Float Acc[3];
} DataIn;

static void gravity_tree_walk_local(const int ipart, Float* acc);

static void interact(const Float, const Float *, const Float, Float *accel);
static bool open_node(const Float, const Float, const Float, const int);

static inline Float node_size(const int node);


void Gravity_Tree_Acceleration()
{
	#pragma omp for
	for (int ipart = 0; ipart < Task.Npart_Total; ipart++) {
			
		Float grav_accel[3] = { 0 };

		gravity_tree_walk_local(ipart, grav_accel);

		P[ipart].Acc[0] += grav_accel[0];
		P[ipart].Acc[1] += grav_accel[1];
		P[ipart].Acc[2] += grav_accel[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
		P[ipart].Grav_Acc[0] = grav_accel[0];
		P[ipart].Grav_Acc[1] = grav_accel[1];
		P[ipart].Grav_Acc[2] = grav_accel[2];
#endif
	} // ipart

	return ;
}
	
static void gravity_tree_walk_local(const int ipart, Float* accel)
{
	int node = 0;

// if (is_local_part)
// copy_part_to_dataout(ipart, D);

	//const Float iAcc = ALENGTH3(P[ipart].Acc);

	while (Tree[node].DNext != 0) {
		
		/* if (is_not_local(node) && is_local_part) {
		 *	
		 *	export_particle(grav_data, Tree[node].DNext);
		 *
		 * 	node++;
		 *
		 *	continue;
		 * }
		 */
	
		Float dr[3] = {P[ipart].Pos[0] - Tree[node].CoM[0],	
			P[ipart].Pos[1] - Tree[node].CoM[1], 
			P[ipart].Pos[2] - Tree[node].CoM[2]};

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		if (Tree[node].DNext < 0) { // encountered particle
			
			interact(nMass, dr, r2, accel);

			node++;
			
			continue;
		}
		
		Float nSize = node_size(node); // now check opening criteria

		if (nSize*nSize > r2 * NODE_OPEN_PARAM) { // BH criterion

			node++;

			continue;
		}

		Float l = 0.6 * nSize;

		if ( (dr[0] < l) || (dr[1] < l) || (dr[2] < l) ) {
		
			node++;

			continue;
		}

		interact(nMass, dr, r2, accel); // use node

		node += fmax(1, Tree[node].DNext);
		
	} // while

	return ;
}

static bool relative_criterion(const Float r2, const Float iAcc, 
		const int node, const Float size)
{
	Float val = Const.Gravity * Tree[node].Mass * size*size / (r2*r2);

	return (val <= NODE_OPEN_PARAM * iAcc);
}

static const double h = GRAV_SOFTENING / 3.0; // Plummer equivalent softening

static void interact(const Float mass, const Float dr[3], const Float r2, 
		Float accel[3])
{
	const Float r = sqrt(r2);
	
	Float rinv = 1/r;

	if (r < h) {
	
		Float u = r/h;
		Float u2 = u*u;
		Float u3 = u2*u;
			
		rinv = sqrt(14*u- 84*u3 + 140*u2*u2 - 90*u2*u3 + 21*u3*u3) / h;
	} 
	
	Float acc_mag = Const.Gravity * mass * p2(rinv);

	accel[0] += -acc_mag * dr[0] * rinv;
	accel[1] += -acc_mag * dr[1] * rinv;
	accel[2] += -acc_mag * dr[2] * rinv;

	return ;
}

static inline Float node_size(const int node)
{
	int lvl = Tree[node].Bitfield & 0x3F; // level

	return Domain.Size / ((Float) (1UL << lvl)); // Domain.Size/2^level
}

#endif // GRAVITY_TREE
