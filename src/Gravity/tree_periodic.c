#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "tree.h"
#include "gravity_periodic.h"

#if defined(GRAVITY_TREE) && defined (PERIODIC)

static struct Walk_Data_Particle copy_send_from(const int ipart);
static void add_recv_to(const int ipart);
static bool interact_with_topnode(const int);
static void interact_with_topnode_particles(const int j);
static void gravity_tree_walk_ewald(const int tree_start);
static void gravity_tree_walk_ewald_BH(const int tree_start);
static void interact_with_ewald_cube(const Float *, const Float);

/*
 * Compute the correction to the gravitational force due the periodic
 * infinite box using the tree and the Ewald method (Hernquist+ 1992).
 * This is widely identical to tree_accel, except for the opening criteria.
 */

static struct Walk_Data_Particle Send = { 0 };
static struct Walk_Data_Result Recv = { 0 };
#pragma omp threadprivate(Send,Recv)

void Gravity_Tree_Periodic()
{
	Profile("Grav Tree Periodic");

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		memset(&Send, 0, sizeof(Send));
		memset(&Recv, 0, sizeof(Recv));

		Send = copy_send_from(ipart);

		for (int j = 0; j < NTop_Nodes; j++) {

				if (interact_with_topnode(j))
					continue;

				//if (D[j].TNode.Target < 0) { // not local ?
				//
				//	export_particle_to_rank(ipart, -target-1);	
				//
				//	continue;
				//}

				if (D[j].TNode.Npart <= VECTOR_SIZE) { // open top leave

					interact_with_topnode_particles(j);

					continue;
				}
			
			int tree_start = D[j].TNode.Target;

			if (Sig.Use_BH_Criterion)  // use BH criterion
				gravity_tree_walk_ewald_BH(tree_start);
			else
				gravity_tree_walk_ewald(tree_start);

		} // for j

		add_recv_to(ipart);

	} // for i

	Profile("Grav Tree Periodic");

	return ;
}

static struct Walk_Data_Particle copy_send_from(const int ipart)
{
	struct Walk_Data_Particle Send = { 0 };

	Send.ID = P.ID[ipart];

	Send.Pos[0] = P.Pos[0][ipart];
	Send.Pos[1] = P.Pos[1][ipart];
	Send.Pos[2] = P.Pos[2][ipart];
	
	Send.Acc = P.Last_Acc_Mag[ipart];

	Send.Mass = P.Mass[ipart];

	return Send;
}

static void add_recv_to(const int ipart)
{
	P.Acc[0][ipart] += Recv.Grav_Acc[0];
	P.Acc[1][ipart] += Recv.Grav_Acc[1];
	P.Acc[2][ipart] += Recv.Grav_Acc[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	P.Grav_Acc[0][ipart] += Recv.Grav_Acc[0];
	P.Grav_Acc[1][ipart] += Recv.Grav_Acc[1];
	P.Grav_Acc[2][ipart] += Recv.Grav_Acc[2];
#endif

#ifdef GRAVITY_POTENTIAL
	P.Grav_Pot[ipart] += Recv.Grav_Potential;
#endif

	P.Cost[ipart] += Recv.Cost;

	return ;
}



static bool interact_with_topnode (const int j)
{
	const Float node_size = Domain.Size / (1UL << D[j].TNode.Level);

	bool want_open_node = false;

	Float dr[3] = { D[j].TNode.CoM[0] - Send.Pos[0],
					D[j].TNode.CoM[1] - Send.Pos[1],
				    D[j].TNode.CoM[2] - Send.Pos[2] };

	Periodic_Nearest(&dr[0]);

	Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

	if (Sig.Use_BH_Criterion) {

		if (node_size*node_size > r2 * TREE_OPEN_PARAM_BH)
			want_open_node = true;

	} else { // relative criterion

		Float node_mass = D[j].TNode.Mass;

		Float fac = Send.Acc/Const.Gravity*TREE_OPEN_PARAM_REL;

		if (node_mass*node_size*node_size > r2*r2 * fac)
			want_open_node = true;
	}

	double dx = fabs(Send.Pos[0] - D[j].TNode.Pos[0]);

	if (dx < 0.6 * node_size) {  // inside subtree ? -> always walk
	
		double dy = fabs(Send.Pos[1] - D[j].TNode.Pos[1]);

		if (dy < 0.6 * node_size) {
		
			double dz = fabs(Send.Pos[2] - D[j].TNode.Pos[2]);
		
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

	interact_with_ewald_cube(dr, D[j].TNode.Mass);

	return true;
}


/*
 * Top nodes with less than 8 particles point not to the tree but to P as 
 * targets. As we have to open this one and it is local, we directly interact
 * with the particles.
 */

static void interact_with_topnode_particles(const int j)
{
	const int first = D[j].TNode.Target;
	const int last = first + D[j].TNode.Npart;

	for (int jpart = first; jpart < last; jpart++) {

		Float dr[3] = { P.Pos[0][jpart] - Send.Pos[0],
					    P.Pos[1][jpart] - Send.Pos[1],
			            P.Pos[2][jpart] - Send.Pos[2] };

		Periodic_Nearest(dr);

		if (dr[0] != 0)
			if(dr[1] != 0)
				if(dr[2] != 0)
					interact_with_ewald_cube(dr, P.Mass[jpart]);
	}

	return ;
}


/*
 * This walks a subtree starting at "tree_start" to estimate the Ewald 
 * correction to the gravitational force. Because the force is proportional to
 * the distance a different opening criterion can be used, similar to 
 * Gadget-2. However, this cannot happen across a periodic boundary. 
 */

static void gravity_tree_walk_ewald(const int tree_start)
{
	const Float fac = Send.Acc / Const.Gravity * TREE_OPEN_PARAM_REL ;

	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		bool want_open_node = false;

		if (Tree[node].DNext < 0) { // particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++) {

				Float dr[3] = { P.Pos[0][jpart] - Send.Pos[0],
							    P.Pos[1][jpart] - Send.Pos[1],
					            P.Pos[2][jpart] - Send.Pos[2] };

				Periodic_Nearest(dr);

				if (dr[0] != 0)
					if (dr[1] != 0) 
						if (dr[2] != 0)
							interact_with_ewald_cube(dr, P.Mass[jpart]);

			} // for jpart

			node++;

			continue;
		}

		Float dr[3] = {Tree[node].CoM[0] - Send.Pos[0],
					    Tree[node].CoM[1] - Send.Pos[1],
					    Tree[node].CoM[2] - Send.Pos[2]};

		Periodic_Nearest(dr);

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float node_mass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		Float ds[3] = { Tree[node].Pos[0]- Send.Pos[0],
						Tree[node].Pos[1]- Send.Pos[1],
						Tree[node].Pos[2]- Send.Pos[2] };	

		if (node_mass*node_size*node_size > r2*r2 * fac) 
			want_open_node = true;
		else if (fabs(ds[0]) < 0.6 * node_size)  // particle inside node ?
				if (fabs(ds[1]) < 0.6 * node_size) 
					if (fabs(ds[2]) < 0.6 * node_size) 
						want_open_node = true;

		Periodic_Nearest(ds);

		if (want_open_node) { // check if we really have to open

			double dl = 0.5 * (Boxsize - node_size);

			if (fabs(ds[0]) > dl) {

				node++;

				continue;
			}

			if (fabs(ds[1]) > dl) {

				node++;

				continue;
			}

			if (fabs(ds[2]) > dl) {

				node++;

				continue;
			}

			if (node_size > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}

		} // if want_open_node

		interact_with_ewald_cube(dr, node_mass);

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}

static void gravity_tree_walk_ewald_BH(const int tree_start)
{
	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (node != tree_end) {

		bool want_open_node = false;

		if (Tree[node].DNext < 0) { // particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++) {

				Float dr[3] = {P.Pos[0][jpart] - Send.Pos[0],
							    P.Pos[1][jpart] - Send.Pos[1],
					            P.Pos[2][jpart] - Send.Pos[2]};

				Periodic_Nearest(dr);

				if (dr[0] != 0)
					if (dr[1] != 0) 
						if (dr[2] != 0)
							interact_with_ewald_cube(dr, P.Mass[jpart]);

			} // for jpart

			node++;

			continue;
		}

		Float dr[3] = {Tree[node].CoM[0] - Send.Pos[0],
					    Tree[node].CoM[1] - Send.Pos[1],
					    Tree[node].CoM[2] - Send.Pos[2]};

		Periodic_Nearest(dr);

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float node_mass = Tree[node].Mass;

		Float node_size = Node_Size(node); // now check opening criteria

		Float ds[3] = { Tree[node].Pos[0]- Send.Pos[0],
						Tree[node].Pos[1]- Send.Pos[1],
						Tree[node].Pos[2]- Send.Pos[2] };	

		if (node_size*node_size > r2 * p2(TREE_OPEN_PARAM_BH)) 
			want_open_node = true;
		else if (fabs(ds[0]) < 0.6 * node_size)  // particle inside node ?
				if (fabs(ds[1]) < 0.6 * node_size) 
					if (fabs(ds[2]) < 0.6 * node_size) 
						want_open_node = true;

		Periodic_Nearest(ds);

		if (want_open_node) { // check if we really have to open

			double dl = 0.5 * (Boxsize - node_size);

			if (fabs(ds[0]) > dl) {

				node++;

				continue;
			}

			if (fabs(ds[1]) > dl) {

				node++;

				continue;
			}

			if (fabs(ds[2]) > dl) {

				node++;

				continue;
			}

			if (node_size > 0.2 * Boxsize) { // too large

				node++;

				continue;
			}

		} // if want_open_node

		interact_with_ewald_cube(dr, node_mass);

		node += Tree[node].DNext; // skip sub-branch, goto sibling

	} // while

	return ;
}




static void interact_with_ewald_cube (const Float * dr, const Float mass)
{
	Float result[3] = { 0 };

	Ewald_Correction(dr, &result[0]);

	Recv.Grav_Acc[0] += Const.Gravity * mass * result[0];
	Recv.Grav_Acc[1] += Const.Gravity * mass * result[1];
	Recv.Grav_Acc[2] += Const.Gravity * mass * result[2];
	
	Recv.Cost++;

#ifdef GRAVITY_POTENTIAL
	
	Float pot_corr = 0;

	Ewald_Potential(dr, &pot_corr);
	
	Recv.Grav_Potential += Const.Gravity * mass * pot_corr;

#endif // GRAVITY_POTENTIAL

	return ;
}


#undef N_EWALD

#endif
