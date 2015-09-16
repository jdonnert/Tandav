#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#ifdef GRAVITY_TREE

struct Walk_Data_Send { // buffer stores input values for tree walk
	int ipart;
	int pos[3];
	int last_acc_mag;
	Float mass;
} Psend = { 0 };
#pragma omp threadprivate(Psend)

struct Walk_Data_Recv { // buffer stores partial results of tree walk
	Float cost;
	Float grav_acc[3];
#ifdef GRAVITY_POTENTIAL
	Float grav_potential;
#endif
} Precv = { 0 };
#pragma omp threadprivate(Precv)

static void prepare_buffers_from(const int);
static void add_Precv_to(const int);
static bool interact_with_topnode(const int);
static void interact_with_topnode_particles(const int);
static void gravity_tree_walk(const int);
static void gravity_tree_walk_first(const int);
static void interact(const Float,const Float *, const Float);

/*
 * Interact with all local topnodes, either with the node directly, or with
 * the particles it contains (small top node). Or, walk the tree and estimate 
 * gravitational acceleration using two different opening criteria. Also open 
 * all nodes containing ipart to avoid large maximum errors. 
 * Barnes & Hut 1984, Springel 2006, Dehnen & Read 2012.
 */

void Gravity_Tree_Acceleration()
{
	Profile("Grav Tree Walk");

	rprintf("Tree acceleration ");

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		memset(&Precv, 0, sizeof(Precv));

		prepare_buffers_from(ipart);

		bool use_BH_criterion = (Psend.last_acc_mag == 0);

		for (int j = 0; j < NTop_Nodes; j++) {

			if (interact_with_topnode(j))
				continue;

			//if (D[j].TNode.Target < 0) { // not local ?
			//
			//	export_particle_to_rank(-target-1);	
			//
			//	continue;
			//}

			if (D[j].TNode.Npart <= 8) { // open top leaf

				interact_with_topnode_particles(j);

				continue;
			}

			int tree_start = D[j].TNode.Target;

			if (use_BH_criterion)
				gravity_tree_walk_first(tree_start);
			else
				gravity_tree_walk(tree_start);
		} // for j

		 add_Precv_to(ipart);

		// work_queue();

	} // for i

	rprintf(" done \n");

	Profile("Grav Tree Walk");

	return ;
}

/*
 * Copy relevant particle variables into the grav_data buffer and 
 */

static void prepare_buffers_from(const int ipart)
{
	memset(&Psend, 0, sizeof(Psend));

	Psend.ipart = ipart;

	Psend.pos[0] = P[ipart].Pos[0];
	Psend.pos[1] = P[ipart].Pos[1];
	Psend.pos[2] = P[ipart].Pos[2];
	
	Psend.last_acc_mag = ALENGTH3(P[ipart].Acc);
	
	Psend.mass = P[ipart].Mass;


	// sink
	
	P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0; // zero sink

	return ;
}

static void add_Precv_to(const int ipart)
{
	P[ipart].Acc[0] += Precv.grav_acc[0];
	P[ipart].Acc[1] += Precv.grav_acc[1];
	P[ipart].Acc[2] += Precv.grav_acc[2];

	P[ipart].Grav_Acc[0] += Precv.grav_acc[0];
	P[ipart].Grav_Acc[1] += Precv.grav_acc[1];
	P[ipart].Grav_Acc[2] += Precv.grav_acc[2];
		
#ifdef GRAVITY_POTENTIAL
	P[ipart].Grav_Pot = Precv.grav_pot;
#endif

	P[ipart].Cost += Precv.cost;

	return ;
}

/*
 * For top nodes far away, we don't have to do a tree walk or send the particle
 * around. Similar to the normal tree walk we first check if the top node 
 * contains the particle and then check the two criteria.
 */

static bool interact_with_topnode(const int j)
{
	const Float node_size = Domain.Size / ((Float)(1UL << D[j].TNode.Level));

	Float dr[3] = {Psend.pos[0] - D[j].TNode.Pos[0],
				   Psend.pos[1] - D[j].TNode.Pos[1],
				   Psend.pos[2] - D[j].TNode.Pos[2]};

	if (fabs(dr[0]) < 0.6 * node_size) // inside subtree ? -> always walk
		if (fabs(dr[1]) < 0.6 * node_size)
			if (fabs(dr[2]) < 0.6 * node_size)
				return false;

	dr[0] = Psend.pos[0] - D[j].TNode.CoM[0];
	dr[1] = Psend.pos[1] - D[j].TNode.CoM[1];
	dr[2] = Psend.pos[2] - D[j].TNode.CoM[2];

	Periodic_Nearest(&dr[0]); // PERIODIC

	Float r2 = ASCALPROD3(dr);

	Float node_mass = D[j].TNode.Mass;

	if (Psend.last_acc_mag == 0) {

		if (p2(node_size) > r2 * TREE_OPEN_PARAM_BH)
			return false;

	} else {

		Float fac = Psend.last_acc_mag/Const.Gravity * TREE_OPEN_PARAM_REL;

		if (node_mass * p2(node_size) > r2*r2 * fac)
			return false;
	}

	interact(node_mass, dr, r2);

	return true;
}

/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets
 */

static void interact_with_topnode_particles(const int j)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		if (jpart == Psend.ipart)
			continue;

		Float dr[3] = { Psend.pos[0] - P[jpart].Pos[0],
						Psend.pos[1] - P[jpart].Pos[1],
						Psend.pos[2] - P[jpart].Pos[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		interact(P[jpart].Mass, dr, r2);
	}

	return ;
}

/*
 * This function walks the local tree and computes the gravitational 
 * acceleration using the relative opening criterion (Springel 2005).
 * If we encounter a particle bundle we interact with all of them.
 */

static void gravity_tree_walk(const int tree_start)
{
	const Float fac = Psend.last_acc_mag/Const.Gravity*TREE_OPEN_PARAM_REL;
	const Float pos_i[3] = {Psend.pos[0], Psend.pos[1], Psend.pos[2]};
	const int ipart = Psend.ipart;

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

				interact(P[jpart].Mass, dr, r2);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float node_mass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		if (node_mass*p2(node_size) > r2*r2 * fac) { // relative criterion

			node++;

			continue;
		}

		Float dx = fabs(pos_i[0] - Tree[node].Pos[0]); // Springel '06: (19)

		if (dx < 0.6 * node_size) {

			Float dy = fabs(pos_i[1] - Tree[node].Pos[1]);

			if (dy < 0.6 * node_size) {

				Float dz = fabs(pos_i[2] - Tree[node].Pos[2]);

				if (dz < 0.6 * node_size) {

					node++;

					continue;
				}
			}
		}

		interact(node_mass, dr, r2); // use node

		node += Tree[node].DNext; // skip branch

	} // while

	return ;
}

/*
 * Walk tree and use the B&H opening criterion, which does not require a prior
 * particle acceleration.
 */

static void gravity_tree_walk_first(const int tree_start)
{
	const Float pos_i[3] = {Psend.pos[0], Psend.pos[1], Psend.pos[2]};

	int node = tree_start;

	while (Tree[node].DNext != 0 || node == tree_start) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -(Tree[node].DNext + 1); // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == Psend.ipart)
					continue;

				Float dr[3] = { pos_i[0] - P[jpart].Pos[0],
							    pos_i[1] - P[jpart].Pos[1],
					            pos_i[2] - P[jpart].Pos[2] };

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				interact(P[jpart].Mass, dr, r2);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float node_mass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		if (node_size*node_size > r2 * TREE_OPEN_PARAM_BH) { // BH criterion

			node++; // open

			continue;
		}

		interact(node_mass, dr, r2); // use node

		node += fmax(1, Tree[node].DNext);

	} // while

	return ;
}


/*
 * Gravitational force law using the K1 softening kernel with central 
 * value corresponding to Plummer softening (Dehnen 2001).
 */

static void interact(const Float mass, const Float dr[3], const Float r2)
{
	//const Float h = 105/32 * GRAV_SOFTENING; // Plummer equiv softening
	const Float h =  GRAV_SOFTENING/3; // Plummer equiv softening

	const Float r = sqrt(r2);

	Float r_inv = 1/r;

#ifdef GRAVITY_POTENTIAL
	Float r_inv_pot = r_inv;
#endif

	if (r < h) { // soften 1/r with Wendland C2 kernel

		Float u = r/h;
		Float u2 = u*u;
		Float u3 = u2*u;

		//r_inv = sqrt(u * (135*u2*u2 - 294*u2 + 175)/(16*h*h) );
		r_inv = sqrt(14*u - 84*u3 + 140*u2*u2 - 90*u2*u3 + 21*u3*u3)/h;
#ifdef GRAVITY_POTENTIAL
		r_inv_pot = ( 45*u3*u3  - 147*u2*u2 + 175*u2 - 105) /(32*h);
#endif
	}

	Float acc_mag = -Const.Gravity * mass * p2(r_inv);

	Precv.grav_acc[0] += acc_mag * dr[0] * r_inv;
	Precv.grav_acc[1] += acc_mag * dr[1] * r_inv;
	Precv.grav_acc[2] += acc_mag * dr[2] * r_inv;

#ifdef GRAVITY_POTENTIAL
	Precv.grav_potential += -Const.Gravity * mass * r_inv_pot;
#endif

	Precv.cost += 1; // count interactions

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
