#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#ifdef GRAVITY_TREE

struct Walk_Data_Send { 
	int Ipart;
	int Pos[3];
	int Acc; 	// magnitude of the last acceleration
	Float Mass;
};

struct Walk_Data_Recv { 
	Float Cost;
	Float Grav_Acc[3];
#ifdef GRAVITY_POTENTIAL
	Float Grav_Potential;
#endif
};

static struct Walk_Data_Send prepare_send_from(const int ipart);
static void write_recv_to(const int ipart, const struct Walk_Data_Recv);

static bool interact_with_topnode(const int, const struct Walk_Data_Send, 
		struct Walk_Data_Recv * restrict);
static void interact_with_topnode_particles(const int, 
		const struct Walk_Data_Send, struct Walk_Data_Recv * restrict);
static void interact(const Float, const Float *, const Float, 
		struct Walk_Data_Recv * restrict);

static void gravity_tree_walk(const int, const struct Walk_Data_Send, 
		struct Walk_Data_Recv * restrict);
static void gravity_tree_walk_first(const int, const struct Walk_Data_Send, 
		struct Walk_Data_Recv * restrict);

/*
 * We do not walk referencing particles, but copy the required particle data
 * into a send buffer "send". The results are written into a sink buffer "recv".
 * This also makes communication easier.
 * These buffers interact with all local topnodes, either with the node 
 * directly, or with the particles it contains (small top node). Or, walk 
 * the tree and estimate gravitational acceleration using two different
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

		struct Walk_Data_Recv recv = { 0 };
		struct Walk_Data_Send send = prepare_send_from(ipart);

		P[ipart].Acc[0] = P[ipart].Acc[1] = P[ipart].Acc[2] = 0;
	
		for (int j = 0; j < NTop_Nodes; j++) {

			if (interact_with_topnode(j, send, &recv))
				continue;

			//if (D[j].TNode.Target < 0) { // not local ?
			//
			//	export_particle_to_rank(ipart, -target-1);	
			//
			//	continue;
			//}

			if (D[j].TNode.Npart <= 8) { // open top leave

				interact_with_topnode_particles(j, send, &recv);

				continue;
			}

			int tree_start = D[j].TNode.Target;

			if (send.Acc == 0) // use BH criterion
				gravity_tree_walk_first(tree_start, send, &recv);
			else
				gravity_tree_walk(tree_start, send, &recv);
		} // for j

		write_recv_to(ipart, recv);

	} // for i

	rprintf(" done \n");

	Profile("Grav Tree Walk");

	#pragma omp barrier

	return ;
}

static struct Walk_Data_Send prepare_send_from(const int ipart)
{
	struct Walk_Data_Send send = { 0 };

	send.Ipart = ipart;

	send.Pos[0] = P[ipart].Pos[0];
	send.Pos[1] = P[ipart].Pos[1];
	send.Pos[2] = P[ipart].Pos[2];
	
	send.Acc = ALENGTH3(P[ipart].Acc);

	send.Mass = P[ipart].Mass;

	return send;
}

static void write_recv_to(const int ipart, const struct Walk_Data_Recv recv)
{
	P[ipart].Acc[0] = recv.Grav_Acc[0];
	P[ipart].Acc[1] = recv.Grav_Acc[1];
	P[ipart].Acc[2] = recv.Grav_Acc[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	P[ipart].Grav_Acc[0] = recv.Grav_Acc[0];
	P[ipart].Grav_Acc[1] = recv.Grav_Acc[1];
	P[ipart].Grav_Acc[2] = recv.Grav_Acc[2];
#endif

#ifdef GRAVITY_POTENTIAL
	P[ipart].Grav_Pot = recv.Grav_Potential;
#endif

	P[ipart].Cost = recv.Cost;

	return ;
}

/*
 * For top nodes far away, we don't have to do a tree walk or send the particle
 * around. Similar to the normal tree walk we first check if the top node 
 * contains the particle and then check the two criteria.
 */

static bool interact_with_topnode(const int j, const struct Walk_Data_Send send,
		struct Walk_Data_Recv * restrict recv)
{
	const Float nSize = Domain.Size / ((Float)(1UL << D[j].TNode.Level));

	Float dr[3] = {send.Pos[0] - D[j].TNode.Pos[0],
				   send.Pos[1] - D[j].TNode.Pos[1],
				   send.Pos[2] - D[j].TNode.Pos[2]};
	
	if (fabs(dr[0]) < 0.6 * nSize) // inside subtree ? -> always walk
		if (fabs(dr[1]) < 0.6 * nSize)
			if (fabs(dr[2]) < 0.6 * nSize)
				return false; 

	dr[0] = send.Pos[0] - D[j].TNode.CoM[0];
	dr[1] = send.Pos[1] - D[j].TNode.CoM[1];
	dr[2] = send.Pos[2] - D[j].TNode.CoM[2];

	Periodic_Nearest(&dr[0]); // PERIODIC

	Float r2 = ASCALPROD3(dr);

	Float node_mass = D[j].TNode.Mass;

	if (send.Acc == 0) {

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH)
			return false;

	} else {

		Float fac = send.Acc/Const.Gravity * TREE_OPEN_PARAM_REL;

		if (node_mass*nSize*nSize > r2*r2 * fac)
			return false;
	}
	
	interact(node_mass, dr, r2, recv);

	return true;
}

/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets
 */

static void interact_with_topnode_particles(const int j, 
		const struct Walk_Data_Send send, 
		struct Walk_Data_Recv * restrict recv)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		if (jpart == send.Ipart)
			continue;

		Float dr[3] = { send.Pos[0] - P[jpart].Pos[0],
					    send.Pos[1] - P[jpart].Pos[1],
			        	send.Pos[2] - P[jpart].Pos[2] };

		Periodic_Nearest(dr); // PERIODIC
		
		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float mpart = P[jpart].Mass;

		interact(mpart, dr, r2, recv);
	}

	return ;
}

/*
 * This function walks the local tree and computes the gravitational 
 * acceleration using the relative opening criterion (Springel 2005).
 * If we encounter a particle bundle we interact with all of them.
 */

static void gravity_tree_walk(const int tree_start, 
		const struct Walk_Data_Send send, struct Walk_Data_Recv * restrict recv)
{
	const Float fac = send.Acc / Const.Gravity * TREE_OPEN_PARAM_REL;

	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == send.Ipart)
					continue;

				Float dr[3] = { send.Pos[0] - P[jpart].Pos[0],
							    send.Pos[1] - P[jpart].Pos[1],
					            send.Pos[2] - P[jpart].Pos[2] };

				Periodic_Nearest(dr); // PERIODIC

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				Float mpart = P[jpart].Mass;

				interact(mpart, dr, r2, recv);
			}

			node++;

			continue;
		}

		Float dr[3] = { send.Pos[0] - Tree[node].CoM[0],
					    send.Pos[1] - Tree[node].CoM[1],
					    send.Pos[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion

			node++;

			continue;
		}

		Float dx = fabs(send.Pos[0] - Tree[node].Pos[0]); // part in node ?

		if (dx < 0.6 * nSize) {  

			Float dy = fabs(send.Pos[1] - Tree[node].Pos[1]);

			if (dy < 0.6 * nSize) {

				Float dz = fabs(send.Pos[2] - Tree[node].Pos[2]);

				if (dz < 0.6 * nSize) {

					node++;

					continue;
				}
			}
		}

		interact(nMass, dr, r2, recv); // use node

		node += Tree[node].DNext; // skip branch

	} // while

	return ;
}

/*
 * Walk tree and use the B&H opening criterion, which does not require a prior
 * particle acceleration.
 */

static void gravity_tree_walk_first(const int tree_start, 
		const struct Walk_Data_Send send, struct Walk_Data_Recv * restrict recv)
{
	int node = tree_start;

	while (Tree[node].DNext != 0 || node == tree_start) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -(Tree[node].DNext + 1); // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == send.Ipart)
					continue;

				Float dr[3] = { send.Pos[0] - P[jpart].Pos[0],
							    send.Pos[1] - P[jpart].Pos[1],
					            send.Pos[2] - P[jpart].Pos[2] };

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				Float mpart = P[jpart].Mass;

				interact(mpart, dr, r2, recv);
			}

			node++;

			continue;
		}

		Float dr[3] = { send.Pos[0] - Tree[node].CoM[0],
					    send.Pos[1] - Tree[node].CoM[1],
					    send.Pos[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH) { // BH criterion

			node++; // open

			continue;
		}

		interact(nMass, dr, r2, recv); // use node

		node += fmax(1, Tree[node].DNext);

	} // while

	return ;
}


/*
 * Gravitational force law using Wendland C2 softening kernel with central 
 * value corresponding to Plummer softening.
 */

static void interact(const Float mass, const Float dr[3], const Float r2, 
		struct Walk_Data_Recv * restrict recv)
{
	//const Float h = GRAV_SOFTENING / 3.0; // Plummer equiv softening
	const Float h = 105/32 * GRAV_SOFTENING;
	const Float r = sqrt(r2);

	Float r_inv = 1/r;

#ifdef GRAVITY_POTENTIAL
	Float r_inv_pot = r_inv;
#endif

	if (r < h) { // soften 1/r with Wendland C2 kernel

		Float u = r/h;
		Float u2 = u*u;
	//	Float u3 = u2*u;

	//	r_inv = sqrt(14*u - 84*u3 + 140*u2*u2 - 90*u2*u3 + 21*u3*u3) / h;
		r_inv = sqrt(u * (135*u2*u2 - 294*u2 + 175))/(4*h) ;
		

#ifdef GRAVITY_POTENTIAL
		r_inv_pot = (7*u2 - 21*u2*u2 + 28*u3*u2 - 15*u3*u3 + u3*u3*u*8 - 3) / h;
#endif
	}

	Float acc_mag = -Const.Gravity * mass * p2(r_inv);

	recv->Grav_Acc[0] += acc_mag * dr[0] * r_inv;
	recv->Grav_Acc[1] += acc_mag * dr[1] * r_inv;
	recv->Grav_Acc[2] += acc_mag * dr[2] * r_inv;

#ifdef GRAVITY_POTENTIAL
	recv->Grav_Potential += -Const.Gravity * mass * r_inv_pot;
#endif

	recv->Cost++;

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
