#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "gravity_tree.h"
#include "gravity_periodic.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
										Float* Accel, Float *Pot);

/*
 * Compute the correction to the gravitational force due the periodic
 * infinite box using the tree and the Ewald method (Hernquist+ 1992).
 */

void Gravity_Tree_Periodic()
{
	Profile("Grav Tree Periodic");

	Profile("Grav Tree Periodic");

	return ;
}

/*
 * This walks a subtree starting at "tree_start" to estimate the Ewald 
 * correction to the gravitational force. Because the force is proportional to
 * the distance a different opening criterion can be used, similar to Gadget-2.
 * However, this cannot happen across a periodic boundary. 
 */

static void gravity_tree_periodic_ewald(const int ipart, const int tree_start,
		Float Accel[3], Float Pot[1])
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

				Ewald_Correction(dr, Accel); 

 				Ewald_Potential(dr, Pot);
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

		Float nSize = Node_Size(node); // now check opening criteria

		if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion
			
			want_open_node = true;

			goto skip;
		}
		
		dr[0] = pos_i[0] - Tree[node].Pos[0];
		dr[1] = pos_i[1] - Tree[node].Pos[1];
		dr[2] = pos_i[2] - Tree[node].Pos[2];

		if (dr[0] < 0.6 * nSize) { // particle inside node ? 
			
			if (dr[1] < 0.6 * nSize) {

				if (dr[2] < 0.6 * nSize) {
				
					want_open_node = true;
				
				}
			}
		}

		skip:;

		if (want_open_node) { // check if we really have to open

			Periodic_Nearest(dr);

			Float size = Node_Size(node);
			
			if (fabs(dr[0]) > 0.5 * (Boxsize - size)) { 
				
				node++;

				continue;
			}

			if (fabs(dr[1]) > 0.5 * (Boxsize - size)) {
				
				node++;

				continue;
			
			}

			if (fabs(dr[2]) > 0.5 * (Boxsize - size)) {
				
				node++;

				continue;
			
			}

			if (Node_Size(node) > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}
		
		} // if want_open_node

		Ewald_Correction(dr, Accel); // use node

 		Ewald_Potential(dr, Pot); // OUTPUT_GRAVITATIONAL_POTENTIAL

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}


/*
 * Do the periodic mapping on a 3D distance array. We rely on link time 
 * optimization of the compiler to do the inlining for us. Make sure to put
 * the appropriate compiler switches.
 */

void Periodic_Nearest(Float dr[3])
{
	for (int i = 0; i < 3; i++) {
	
		if (dr[i] > Boxhalf)
			dr[i] -= Boxsize;
		else if (dr[i] < -Boxhalf)
			dr[i] += Boxsize;
	}

	return ;
}

#undef N_EWALD

#endif
