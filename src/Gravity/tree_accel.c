#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#ifdef GRAVITY_TREE

static bool interact_with_topnode(const int, const int, Float*, Float*);
static void interact_with_topnode_particles(const int, const int,
											Float*, Float*);
static void gravity_tree_walk(const int , const int, Float*, Float*);
static void gravity_tree_walk_first(const int , const int, Float*, Float*);
static void interact(const Float,const Float *,const Float,Float *,Float *);

/*
 * Interact with all local topnodes, either with the node directly, or with the 
 * particles it contains (small top node). Or, walk the tree and estimate 
 * gravitational acceleration using two different
 * opening criteria. Also open all nodes containing ipart to avoid large 
 * maximum errors. Barnes & Hut 1984, Springel 2006, Dehnen & Read 2012.
 */

void Gravity_Tree_Acceleration()
{
	Profile("Grav Tree Walk");

	rprintf("Tree acceleration ");

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		bool use_BH_criterion = (ASCALPROD3(P[ipart].Acc) == 0);

		Float grav_accel[3] = { 0 };
		Float pot = 0;

		for (int j = 0; j < NTop_Nodes; j++) {

			if (interact_with_topnode(ipart, j, &grav_accel[0], &pot))
				continue;

			//if (D[j].TNode.Target < 0) { // not local ?
			//
			//	export_particle_to_rank(ipart, -target-1);	
			//
			//	continue;
			//}

			if (D[j].TNode.Npart <= 8) { // open top leave

				interact_with_topnode_particles(ipart, j, grav_accel, &pot);

				continue;
			}

			int tree_start = D[j].TNode.Target;

			if (use_BH_criterion)
				gravity_tree_walk_first(ipart, tree_start, grav_accel, &pot);
			else
				gravity_tree_walk(ipart, tree_start, grav_accel, &pot);
		} // for j

		P[ipart].Acc[0] = grav_accel[0];
		P[ipart].Acc[1] = grav_accel[1];
		P[ipart].Acc[2] = grav_accel[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
		P[ipart].Grav_Acc[0] = grav_accel[0];
		P[ipart].Grav_Acc[1] = grav_accel[1];
		P[ipart].Grav_Acc[2] = grav_accel[2];
#endif

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = pot;
#endif

	} // for i

	rprintf(" done \n");

	Profile("Grav Tree Walk");

	#pragma omp barrier

	return ;
}

/*
 * For top nodes far away, we don't have to do a tree walk or send the particle
 * around. Similar to the normal tree walk we first check if the top node 
 * contains the particle and then check the two criteria.
 */

static bool interact_with_topnode(const int ipart, const int j,
									Float *grav_accel, Float *pot)
{
	const Float nSize = Domain.Size / ((Float)(1UL << D[j].TNode.Level));

	Float dr[3] = {P[ipart].Pos[0] - D[j].TNode.Pos[0],
				   P[ipart].Pos[1] - D[j].TNode.Pos[1],
				   P[ipart].Pos[2] - D[j].TNode.Pos[2]};
	
	if (fabs(dr[0]) < 0.6 * nSize) // inside subtree ? -> always walk
		if (fabs(dr[1]) < 0.6 * nSize)
			if (fabs(dr[2]) < 0.6 * nSize)
				return false; 

	dr[0] = P[ipart].Pos[0] - D[j].TNode.CoM[0];
	dr[1] = P[ipart].Pos[1] - D[j].TNode.CoM[1];
	dr[2] = P[ipart].Pos[2] - D[j].TNode.CoM[2];

	Periodic_Nearest(&dr[0]); // PERIODIC

	Float r2 = ASCALPROD3(dr);

	Float node_mass = D[j].TNode.Mass;

	if (ASCALPROD3(P[ipart].Acc) == 0) {

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH)
			return false;

	} else {

		Float fac = ALENGTH3(P[ipart].Acc)/Const.Gravity * TREE_OPEN_PARAM_REL;

		if (node_mass*nSize*nSize > r2*r2 * fac)
			return false;
	}
	
	interact(node_mass, dr, r2, grav_accel, pot);

	return true;
}

/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets
 */

static void interact_with_topnode_particles(const int ipart, const int j,
		Float *grav_accel, Float *pot)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		if (jpart == ipart)
			continue;

		Float dr[3] = { P[ipart].Pos[0] - P[jpart].Pos[0],
					    P[ipart].Pos[1] - P[jpart].Pos[1],
			            P[ipart].Pos[2] - P[jpart].Pos[2] };

		Periodic_Nearest(dr); // PERIODIC
		
		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float mpart = P[jpart].Mass;

		interact(mpart, dr, r2, grav_accel, pot);
	}

	return ;
}

/*
 * This function walks the local tree and computes the gravitational 
 * acceleration using the relative opening criterion (Springel 2005).
 * If we encounter a particle bundle we interact with all of them.
 */

static void gravity_tree_walk(const int ipart, const int tree_start,
		Float Accel[3], Float Pot[1])
{
	const Float fac = ALENGTH3(P[ipart].Acc) / Const.Gravity
													* TREE_OPEN_PARAM_REL;

	const Float pos_i[3] = {P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2]};

	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == ipart)
					continue;

				Float dr[3] = { pos_i[0] - P[jpart].Pos[0],
							    pos_i[1] - P[jpart].Pos[1],
					            pos_i[2] - P[jpart].Pos[2] };

				Periodic_Nearest(dr); // PERIODIC

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				Float mpart = P[jpart].Mass;

				interact(mpart, dr, r2, Accel, Pot);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion

			node++;

			continue;
		}

		Float dx = fabs(pos_i[0] - Tree[node].Pos[0]); // particle in node ?

		if (dx < 0.6 * nSize) {  

			Float dy = fabs(pos_i[1] - Tree[node].Pos[1]);

			if (dy < 0.6 * nSize) {

				Float dz = fabs(pos_i[2] - Tree[node].Pos[2]);

				if (dz < 0.6 * nSize) {

					node++;

					continue;
				}
			}
		}

		interact(nMass, dr, r2, Accel, Pot); // use node

		node += Tree[node].DNext; // skip branch

	} // while

	return ;
}

/*
 * Walk tree and use the B&H opening criterion, which does not require a prior
 * particle acceleration.
 */

static void gravity_tree_walk_first(const int ipart, const int tree_start,
		Float* Accel, Float *Pot)
{
	const Float pos_i[3] = {P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2]};

	int node = tree_start;

	while (Tree[node].DNext != 0 || node == tree_start) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -(Tree[node].DNext + 1); // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == ipart)
					continue;

				Float dr[3] = { pos_i[0] - P[jpart].Pos[0],
							    pos_i[1] - P[jpart].Pos[1],
					            pos_i[2] - P[jpart].Pos[2] };

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				Float mpart = P[jpart].Mass;

				interact(mpart, dr, r2, Accel, Pot);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH) { // BH criterion

			node++; // open

			continue;
		}

		interact(nMass, dr, r2, Accel, Pot); // use node

		node += fmax(1, Tree[node].DNext);

	} // while

	return ;
}


/*
 * Gravitational force law using Wendland C2 softening kernel with central 
 * value corresponding to Plummer softening.
 */

static void interact(const Float mass, const Float dr[3], const Float r2,
		Float Accel[3], Float *Pot)
{
	const Float h = GRAV_SOFTENING / 3.0; // Plummer equiv softening

	const Float r = sqrt(r2);

	Float r_inv = 1/r;

#ifdef GRAVITY_POTENTIAL
	Float r_inv_pot = r_inv;
#endif

	if (r < h) { // soften 1/r with Wendland C2 kernel

		Float u = r/h;
		Float u2 = u*u;
		Float u3 = u2*u;

		r_inv = sqrt(14*u - 84*u3 + 140*u2*u2 - 90*u2*u3 + 21*u3*u3) / h;

#ifdef GRAVITY_POTENTIAL
		r_inv_pot = (7*u2 - 21*u2*u2 + 28*u3*u2 - 15*u3*u3 + u3*u3*u*8 - 3) / h;
#endif
	}

	Float acc_mag = -Const.Gravity * mass * p2(r_inv);

	Accel[0] += acc_mag * dr[0] * r_inv;
	Accel[1] += acc_mag * dr[1] * r_inv;
	Accel[2] += acc_mag * dr[2] * r_inv;

#ifdef GRAVITY_POTENTIAL
	*Pot += -Const.Gravity * mass * r_inv_pot;
#endif

	return ;
}

/*
 * Bitfield functions on global Tree
 */

Float Node_Size(const int node)
{
	int lvl = Tree[node].Bitfield & 0x3F; // level

	return Domain.Size / ((Float) (1ULL << lvl)); // Domain.Size/2^level
}

int Level(const int node)
{
	return Tree[node].Bitfield & 0x3FUL; // return bit 0-5
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

#endif // GRAVITY_TREE
