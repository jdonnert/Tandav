#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
										Float* Accel, Float *Pot);
static bool interact_with_topnode(const int ipart, const int j,
									Float *fr, Float *fp);
static void interact_with_topnode_particles(const int ipart, const int j,
		Float *fr, Float *fp);

/*
 * Compute the correction to the gravitational force due the periodic
 * infinite box using the tree and the Ewald method (Hernquist+ 1992).
 * This is widely identical to tree_accel, except for the opening criteria.
 */

void Gravity_Tree_Periodic()
{
	Profile("Grav Tree Periodic");

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];
		
		Float ewald_corr[3] = { 0 };
		Float ewald_pot = 0;

		for (int j = 0; j < NTop_Nodes; j++) {
					
			Float fr[3] = { 0 };
			Float fp = 0;

			if (interact_with_topnode(ipart, j, &fr[0], &fp))
				continue;

			//if (D[j].TNode.Target < 0) {
			//
			//	export_particle_to_rank(ipart, -target-1);	
			//
			//	continue;
			//}
			
			if (D[j].TNode.Npart <= 8) { // open top leave

				interact_with_topnode_particles(ipart, j, &fr[0], &fp);
		
				continue;
			}

			int tree_start = D[j].TNode.Target;

			gravity_tree_periodic_ewald(ipart, tree_start, &fr[0], &fp);

			ewald_corr[0] += fr[0];			
			ewald_corr[1] += fr[1];			
			ewald_corr[2] += fr[2];			

			ewald_pot += fp;
		} // for j

		P[ipart].Acc[0] += ewald_corr[0];
		P[ipart].Acc[1] += ewald_corr[1];
		P[ipart].Acc[2] += ewald_corr[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
		P[ipart].Grav_Acc[0] += ewald_corr[0];
		P[ipart].Grav_Acc[1] += ewald_corr[1];
		P[ipart].Grav_Acc[2] += ewald_corr[2];
#endif

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot += ewald_pot;
#endif
	} // for i


	Profile("Grav Tree Periodic");

	return ;
}

static bool interact_with_topnode(const int ipart, const int j,
								 Float *fr, Float *fp)
{
	const Float node_size = Domain.Size / ((Float)(1UL << D[j].TNode.Level));

	bool want_open_node = false;

	Float dr[3] = {P[ipart].Pos[0] - D[j].TNode.CoM[0],
				   P[ipart].Pos[1] - D[j].TNode.CoM[1],
				   P[ipart].Pos[2] - D[j].TNode.CoM[2]};

	Periodic_Nearest(&dr[0]);

	Float r2 = ASCALPROD3(dr);

	if (ASCALPROD3(P[ipart].Acc) == 0) {

		if (node_size*node_size > r2 * TREE_OPEN_PARAM_BH)
			want_open_node = true;

	} else {

		Float nMass = D[j].TNode.Mass;

		Float fac = ALENGTH3(P[ipart].Acc)/Const.Gravity * TREE_OPEN_PARAM_REL;

		if (nMass*node_size*node_size > r2*r2 * fac)
			want_open_node = true;
	}
	
	dr[0] = P[ipart].Pos[0] - D[j].TNode.Pos[0];
	dr[1] = P[ipart].Pos[1] - D[j].TNode.Pos[1];
	dr[2] = P[ipart].Pos[2] - D[j].TNode.Pos[2];

	if (fabs(dr[0]) < 0.6 * node_size) // inside subtree ? -> always walk
		if (fabs(dr[1]) < 0.6 * node_size)
			if (fabs(dr[2]) < 0.6 * node_size) 
				want_open_node = true; 

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

	dr[0] = P[ipart].Pos[0] - D[j].TNode.CoM[0];
	dr[1] = P[ipart].Pos[1] - D[j].TNode.CoM[1];
	dr[2] = P[ipart].Pos[2] - D[j].TNode.CoM[2];

	Periodic_Nearest(dr);

	Ewald_Correction(dr, &fr[0]); // use node

 	Ewald_Potential(dr, &fp[0]); // OUTPUT_GRAVITATIONAL_POTENTIAL

	return true;
}

/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets
 */

static void interact_with_topnode_particles(const int ipart, const int j,
		Float *fr, Float *fp)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		if (jpart == ipart)
			continue;

		Float dr[3] = { P[ipart].Pos[0] - P[jpart].Pos[0],
					    P[ipart].Pos[1] - P[jpart].Pos[1],
			            P[ipart].Pos[2] - P[jpart].Pos[2] };

		Periodic_Nearest(dr);
		
		Ewald_Correction(dr, &fr[0]); // use node

 		Ewald_Potential(dr, &fp[0]); // OUTPUT_GRAVITATIONAL_POTENTIAL
	}

	return ;
}


/*
 * This walks a subtree starting at "tree_start" to estimate the Ewald 
 * correction to the gravitational force. Because the force is proportional to
 * the distance a different opening criterion can be used, similar to Gadget-2.
 * However, this cannot happen across a periodic boundary. 
 */

static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
		Float fr[3], Float fp[1])
{
	const Float fac = ALENGTH3(P[ipart].Acc) / Const.Gravity
		* TREE_OPEN_PARAM_REL;
	
	const Float pos_i[3] = {P[ipart].Pos[0], P[ipart].Pos[1], P[ipart].Pos[2]};
	
	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		bool want_open_node = false;

		if (Tree[node].DNext < 0) { // particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++) {

				if (jpart == ipart)
					continue;

				Float dr[3] = { pos_i[0] - P[jpart].Pos[0],
							    pos_i[1] - P[jpart].Pos[1],
					            pos_i[2] - P[jpart].Pos[2] };

				Periodic_Nearest(dr);

				Ewald_Correction(dr, &fr[0]); 

 				Ewald_Potential(dr, &fp[0]);
			}

			node++;

			continue;
		}

		Float dr[3] = { pos_i[0] - Tree[node].CoM[0],
					    pos_i[1] - Tree[node].CoM[1],
					    pos_i[2] - Tree[node].CoM[2] };

		Periodic_Nearest(dr);

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		if (nMass*node_size*node_size > r2*r2 * fac) // relative criterion
			want_open_node = true;

		dr[0] = pos_i[0] - Tree[node].Pos[0];
		dr[1] = pos_i[1] - Tree[node].Pos[1];
		dr[2] = pos_i[2] - Tree[node].Pos[2];

		if (fabs(dr[0]) < 0.6 * node_size) { // particle inside node ? 
			
			if (fabs(dr[1]) < 0.6 * node_size) {

				if (fabs(dr[2]) < 0.6 * node_size) {
				
					want_open_node = true;
				
				}
			}
		}

		if (want_open_node) { // check if we really have to open

			if (fabs(dr[0]) > 0.5 * (Boxsize - node_size)) { 
				
				node++;

				continue;
			}

			if (fabs(dr[1]) > 0.5 * (Boxsize - node_size)) {
				
				node++;

				continue;
			}

			if (fabs(dr[2]) > 0.5 * (Boxsize - node_size)) {
				
				node++;

				continue;
			}

			if (node_size > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}
		
		} // if want_open_node

		Ewald_Correction(dr, &fr[0]); // use node

 		Ewald_Potential(dr, &fp[0]); // OUTPUT_GRAVITATIONAL_POTENTIAL

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}

#undef N_EWALD

#endif
