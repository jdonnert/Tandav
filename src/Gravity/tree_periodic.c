#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

static bool interact_with_topnode(const int, const struct Walk_Data_Particle, 
		struct Walk_Data_Result * restrict);

static void interact_with_topnode_particles(const int j,
		const struct Walk_Data_Particle, struct Walk_Data_Result * restrict);

static void gravity_tree_walk_ewald(const int tree_start,
		const struct Walk_Data_Particle, struct Walk_Data_Result * restrict);

static void interact_with_ewald_cube(const double *, const Float,
		struct Walk_Data_Result * restrict);

/*
 * Compute the correction to the gravitational force due the periodic
 * infinite box using the tree and the Ewald method (Hernquist+ 1992).
 * This is widely identical to tree_accel, except for the opening criteria.
 */

void Gravity_Tree_Periodic(const struct Walk_Data_Particle send, 
		struct Walk_Data_Result *recv)
{
	for (int j = 0; j < NTop_Nodes; j++) {

			//if (interact_with_topnode(j, send, &recv))
		//		continue;

		//	if (D[j].TNode.Npart <= 8) { // open top leave
//
//				interact_with_topnode_particles(j, send, &recv);
//
//				continue;
//			}
			
		int tree_start = D[j].TNode.Target;

		gravity_tree_walk_ewald(tree_start, send, recv);

	} // for j
	
	return ;
}

static bool interact_with_topnode (const int j,  
		const struct Walk_Data_Particle send, 
		struct Walk_Data_Result * restrict recv )
{
	const Float node_size = Domain.Size / (1UL << D[j].TNode.Level);

	bool want_open_node = false;

	double dr[3] = {D[j].TNode.CoM[0] - send.Pos[0],
					D[j].TNode.CoM[1] - send.Pos[1],
				    D[j].TNode.CoM[2] - send.Pos[2]};

	Periodic_Nearest(&dr[0]);

	Float r2 = ASCALPROD3(dr);

	if (Sig.Use_BH_Criterion) {

		if (node_size*node_size > r2 * TREE_OPEN_PARAM_BH)
			want_open_node = true;

	} else { // relative criterion

		Float nMass = D[j].TNode.Mass;

		Float fac = send.Acc/Const.Gravity*TREE_OPEN_PARAM_REL;

		if (nMass*node_size*node_size > r2*r2 * fac)
			want_open_node = true;
	}

	double dx = fabs(send.Pos[0] - D[j].TNode.Pos[0]);

	if (dx < 0.6 * node_size) {  // inside subtree ? -> always walk
	
		double dy = fabs(send.Pos[1] - D[j].TNode.Pos[1]);

		if (dy < 0.6 * node_size) {
		
			double dz = fabs(send.Pos[2] - D[j].TNode.Pos[2]);
		
			if (dz < 0.6 * node_size) {
			
				want_open_node = true;
			}
		}
	}

	if (want_open_node) { // check if we really have to open

			if (fabs(dr[0]) > 0.5 * (Boxsize - node_size))
				return false;

			if (fabs(dr[1]) > 0.5 * (Boxsize - node_size))
				return false;

			if (fabs(dr[2]) > 0.5 * (Boxsize - node_size))
				return false;

			if (node_size > 0.2 * Boxsize)
				return false;

	} // if want_open_node

	interact_with_ewald_cube(dr, D[j].TNode.Mass, recv);

	return true;
}

/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets. As we have to open this one and it is local, we directly interact
 * with the particles.
 */

static void interact_with_topnode_particles(const int j,
		const struct Walk_Data_Particle send,
		struct Walk_Data_Result * restrict recv)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		double dr[3] = {P[jpart].Pos[0] - send.Pos[0],
					    P[jpart].Pos[1] - send.Pos[1],
			            P[jpart].Pos[2] - send.Pos[2] };

		Periodic_Nearest(dr);
	
		interact_with_ewald_cube(dr, P[jpart].Mass, recv);
	}

	return ;
}


/*
 * This walks a subtree starting at "tree_start" to estimate the Ewald 
 * correction to the gravitational force. Because the force is proportional to
 * the distance a different opening criterion can be used, similar to 
 * Gadget-2. However, this cannot happen across a periodic boundary. 
 */

static void gravity_tree_walk_ewald(const int tree_start,
		const struct Walk_Data_Particle send,
		struct Walk_Data_Result * restrict recv)
{
	const Float fac = send.Acc / (Const.Gravity) * TREE_OPEN_PARAM_REL;

	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		bool want_open_node = false;

		if (Tree[node].DNext < 0) { // particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++) {

				double dr[3] = {P[jpart].Pos[0] - send.Pos[0],
							    P[jpart].Pos[1] - send.Pos[1],
					            P[jpart].Pos[2] - send.Pos[2]};

				Periodic_Nearest(dr);

//printf("%d %g %g %g %g %g %g ", jpart, 
//		P[jpart].Pos[0], P[jpart].Pos[1], P[jpart].Pos[2], 
//	send.Pos[0],send.Pos[1],send.Pos[2]);

				interact_with_ewald_cube(dr, P[jpart].Mass, recv);
			} // for jpart

			node++;

			continue;
		}

		double dr[3] = {Tree[node].CoM[0] - send.Pos[0],
					    Tree[node].CoM[1] - send.Pos[1],
					    Tree[node].CoM[2] - send.Pos[2]};

		Periodic_Nearest(dr);

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		if (Sig.Use_BH_Criterion) { // use BH criterion
	
			if (node_size*node_size > r2 * p2(TREE_OPEN_PARAM_BH)) 
				want_open_node = true;

		} else { // relative criterion

			if (nMass*node_size*node_size > r2*r2 * fac) 
				want_open_node = true;
		}

		double ds[3] = {Tree[node].Pos[0]- send.Pos[0],
						Tree[node].Pos[1]- send.Pos[1],
						Tree[node].Pos[2]- send.Pos[2]};	

		if (fabs(ds[0]) < 0.6 * node_size)  // particle inside node ?
			if (fabs(ds[1]) < 0.6 * node_size) 
				if (fabs(ds[2]) < 0.6 * node_size) 
					want_open_node = true;

		Periodic_Nearest(ds);

		if (want_open_node) { // check if we really have to open

			if (fabs(ds[0]) > 0.5 * (Boxsize - node_size)) {

				node++;

				continue;
			}

			if (fabs(ds[1]) > 0.5 * (Boxsize - node_size)) {

				node++;

				continue;
			}

			if (fabs(ds[2]) > 0.5 * (Boxsize - node_size)) {

				node++;

				continue;
			}

			if (node_size > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}

		} // if want_open_node

		interact_with_ewald_cube(dr, nMass, recv);

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}


static void interact_with_ewald_cube (const double * dr, const Float mass,
		struct Walk_Data_Result * restrict recv)
{
	recv->Cost++;

	Float result[3] = { 0 };

	Ewald_Correction(dr, &result[0]);

	recv->Grav_Acc[0] += Const.Gravity * mass * result[0];
	recv->Grav_Acc[1] += Const.Gravity * mass * result[1];
	recv->Grav_Acc[2] += Const.Gravity * mass * result[2];
	
//printf("%g %g %g %g %g %g %g \n",mass, dr[0], dr[1], dr[2], result[0], result[1], result[2]);
#ifdef GRAVITY_POTENTIAL
	
	Float pot_corr = 0;

	Ewald_Potential(dr, &pot_corr);
	
	recv->Grav_Potential += Const.Gravity * mass * pot_corr;

#endif // GRAVITY_POTENTIAL

	return ;
}


#undef N_EWALD

#endif
