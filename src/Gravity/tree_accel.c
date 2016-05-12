#include "../globals.h"
#include "../domain.h"
#include "gravity.h"
#include "tree.h"
#include "gravity_periodic.h"

#ifdef GRAVITY_TREE

static struct Walk_Data_Particle copy_send_from(const int ipart);
static void add_recv_to(const int ipart);

static bool interact_with_topnode(const int);
static void interact_with_topnode_particles(const int);
static void interact(const Float, const Float *, const Float);

static void gravity_tree_walk(const int);
static void gravity_tree_walk_BH(const int);

static void check_total_momentum(const bool show_change);

/*
 * We do not walk referencing particles, but copy the required particle data
 * into a Send buffer "Send". The results are written into a sink buffer 
 * "Recv". 
 * These buffers interact with all local topnodes, either with the node 
 * directly, or with the particles it contains (small top node). Or, walk 
 * the tree and estimate gravitational acceleration using two different
 * opening criteria. Also open all nodes containing ipart to avoid large 
 * maximum errors. Barnes & Hut 1984, Springel 2006, Dehnen & Read 2012.
 */

static struct Walk_Data_Particle Send = { 0 };
static struct Walk_Data_Result Recv = { 0 };
#pragma omp threadprivate(Send,Recv)

void Gravity_Tree_Acceleration()
{
	Profile("Grav Tree Walk");

	rprintf("Tree acceleration ");

	check_total_momentum(false);

	#pragma omp for schedule(dynamic)
	for (int i = 0; i < NActive_Particles; i++) {

		int ipart = Active_Particle_List[i];

		memset(&Send, 0, sizeof(Send));
		memset(&Recv, 0, sizeof(Recv));

		Send = copy_send_from(ipart);

		P.Last_Acc_Mag[ipart] = sqrt( p2(P.Acc[0][ipart]) 
									+ p2(P.Acc[1][ipart]) 
									+ p2(P.Acc[2][ipart]));

		P.Acc[0][ipart] = P.Acc[1][ipart] = P.Acc[2][ipart] = 0;

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

			if (Sig.Use_BH_Criterion) 
				gravity_tree_walk_BH(tree_start);
			else
				gravity_tree_walk(tree_start);

		} // for j
	
		add_recv_to(ipart);

	} // for i

	Profile("Grav Tree Walk");

	Gravity_Tree_Periodic(); // PERIODIC , add Ewald correction

	rprintf(" done \n");

	check_total_momentum(true);

	return ;
}

static struct Walk_Data_Particle copy_send_from(const int ipart)
{
	struct Walk_Data_Particle Send = { 0 };

	Send.ID = P.ID[ipart];

	Send.Pos[0] = P.Pos[0][ipart];
	Send.Pos[1] = P.Pos[1][ipart];
	Send.Pos[2] = P.Pos[2][ipart];
	
	Send.Acc = sqrt(p2(P.Acc[0][ipart]) + p2(P.Acc[1][ipart]) 
				   + p2(P.Acc[2][ipart]));

	Send.Mass = P.Mass[ipart];

	return Send;
}

static void add_recv_to(const int ipart)
{
	P.Acc[0][ipart] += Recv.Grav_Acc[0];
	P.Acc[1][ipart] += Recv.Grav_Acc[1];
	P.Acc[2][ipart] += Recv.Grav_Acc[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
	P.Grav_Acc[0][ipart] = Recv.Grav_Acc[0];
	P.Grav_Acc[1][ipart] = Recv.Grav_Acc[1];
	P.Grav_Acc[2][ipart] = Recv.Grav_Acc[2];
#endif

#ifdef GRAVITY_POTENTIAL
	P.Grav_Pot[ipart] = Recv.Grav_Potential;
#endif

	P.Cost[ipart] = Recv.Cost;

	return ;
}

/*
 * For top nodes far away, we don't have to do a tree walk or Send the particle
 * around. Similar to the normal tree walk we first check if the top node 
 * contains the particle and then check the two criteria.
 */

static bool interact_with_topnode(const int j)
{
	const Float nSize = Domain.Size / ((Float)(1UL << D[j].TNode.Level));

	Float dr[3] = {D[j].TNode.Pos[0] - Send.Pos[0],
				   D[j].TNode.Pos[1] - Send.Pos[1],
				   D[j].TNode.Pos[2] - Send.Pos[2]};
	
	if (fabs(dr[0]) < 0.6 * nSize) // inside subtree ? -> always walk
		if (fabs(dr[1]) < 0.6 * nSize)
			if (fabs(dr[2]) < 0.6 * nSize)
				return false; 

	dr[0] = D[j].TNode.CoM[0] - Send.Pos[0];
	dr[1] = D[j].TNode.CoM[1] - Send.Pos[1];
	dr[2] = D[j].TNode.CoM[2] - Send.Pos[2];

	Periodic_Nearest(dr); // PERIODIC

	Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

	Float node_mass = D[j].TNode.Mass;

	if (Sig.Use_BH_Criterion) {

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH)
			return false;

	} else {

		Float fac = Send.Acc/Const.Gravity * TREE_OPEN_PARAM_REL;

		if (node_mass*nSize*nSize > r2*r2 * fac)
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

		Float dr[3] = {P.Pos[0][jpart] - Send.Pos[0],
					   P.Pos[1][jpart] - Send.Pos[1] ,
			           P.Pos[2][jpart] - Send.Pos[2] };

		Periodic_Nearest(dr); // PERIODIC
		
		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		if (r2 != 0)
			interact(P.Mass[jpart], dr, r2);
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
	const Float fac = Send.Acc / Const.Gravity * TREE_OPEN_PARAM_REL;

	const int tree_end = tree_start + Tree[tree_start].DNext;

	int node = tree_start;

	while (Tree[node].DNext != 0 || node == tree_start) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				Float dr[3] = {P.Pos[0][jpart] - Send.Pos[0],
								P.Pos[1][jpart] - Send.Pos[1],
								P.Pos[2][jpart] - Send.Pos[2]};
				
				Periodic_Nearest(dr); // PERIODIC

				Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
				
				if (r2 != 0) 
					interact(P.Mass[jpart], dr, r2);
			}

			node++;

			continue;
		}

		Float dr[3] = {Tree[node].CoM[0] - Send.Pos[0],
					   Tree[node].CoM[1] - Send.Pos[1],
					   Tree[node].CoM[2] - Send.Pos[2]};
		
		Periodic_Nearest(dr); // PERIODIC

		Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion
			
			node++;

			continue;
		}

		Float ds[3] = {Tree[node].Pos[0] - Send.Pos[0],
					   Tree[node].Pos[1] - Send.Pos[1],
					   Tree[node].Pos[2] - Send.Pos[2]};

		Periodic_Nearest(ds); // PERIODIC

		if (fabs(ds[0]) < 0.6 * nSize) {  

			if (fabs(ds[1]) < 0.6 * nSize) {

				if (fabs(ds[2]) < 0.6 * nSize) {

					node++;

					continue;
				}
			}
		}

		interact(nMass, dr, r2); // use node

		node += Tree[node].DNext; // skip branch

	} // while

	return ;
}

/*
 * Walk tree and use the B&H opening criterion, which does not require a prior
 * particle acceleration.
 */

static void gravity_tree_walk_BH(const int tree_start)
{
	int node = tree_start;

	while (Tree[node].DNext != 0 || node == tree_start) {

		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -Tree[node].DNext - 1; // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				Float dr[3] = { P.Pos[0][jpart] - Send.Pos[0],
							    P.Pos[1][jpart] - Send.Pos[1],
					            P.Pos[2][jpart] - Send.Pos[2] };

				Periodic_Nearest(dr); // PERIODIC

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				if (r2 != 0)
					interact(P.Mass[jpart], dr, r2);
			}

			node++;

			continue;
		}

		Float dr[3] = {Tree[node].CoM[0] - Send.Pos[0],
					   Tree[node].CoM[1] - Send.Pos[1],
					   Tree[node].CoM[2] - Send.Pos[2]};

		Periodic_Nearest(dr); // PERIODIC

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = Node_Size(node); // now check opening criteria

		if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH) { // BH criterion

			node++; // open

			continue;
		}

		interact(nMass, dr, r2); // use node

		node += Tree[node].DNext;

	} // while

	return ;
}


/*
 * Gravitational force law using Dehnens K1 softening kernel with central 
 * value corresponding to Plummer softening of potential : 
 * h_K1 = -41.0/32.0 * eps_plummer;
 */

static void interact(const Float mass, const Float dr[3], const Float r2)
{
	Float fac = Const.Gravity * mass;
	Float fac_pot = Const.Gravity * mass;

	if (r2 < Epsilon2[1]) { 

		Float u2 = r2 / Epsilon2[1];

		fac *= (175 - u2 * (294 - u2 * 135)) / (16*Epsilon3[1]) ;

		fac_pot *= (u2 * (175 - (u2 * 147  - u2 * 45)) - 105)/(32*Epsilon[1]);

	} else {

		Float r_inv = 1/SQRT(r2); // tempt the compiler to use rsqrtss

		fac *= r_inv * r_inv * r_inv;
		fac_pot *= r_inv;
	}

	Recv.Grav_Acc[0] += fac * dr[0];
	Recv.Grav_Acc[1] += fac * dr[1];
	Recv.Grav_Acc[2] += fac * dr[2];

#ifdef GRAVITY_POTENTIAL
	Recv.Grav_Potential += fac_pot;
#endif

	Recv.Cost += 1;

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

/*
 * Here start MPI communication variables and routines. 
 */


/*
 * Compute total momentum to check the gravity interaction. 
 */

static double Px = 0, Py = 0, Pz = 0, Last = 0;

static void check_total_momentum(const bool show_change)
{
	const int last_DM_part = Task.Npart[0]+Task.Npart[1];

	#pragma omp for reduction(+: Px, Py, Pz)
	for (int ipart = Task.Npart[0]; ipart < last_DM_part; ipart++) {
	
		Px += P.Mass[ipart] * P.Vel[0][ipart];
		Py += P.Mass[ipart] * P.Vel[1][ipart];
		Pz += P.Mass[ipart] * P.Vel[2][ipart];
	}
	
	double ptotal = sqrt( Px*Px + Py*Py + Pz*Pz );

	#pragma omp single
	MPI_Reduce(MPI_IN_PLACE, &ptotal, 1, MPI_DOUBLE, MPI_SUM, MASTER, 
			MPI_COMM_WORLD);

	double rel_err = (ptotal - Last) / Last;

	if (show_change)
		rprintf("Total change in momentum due to gravity: %g \n", rel_err);

	#pragma omp single
	Last = ptotal;

	return ;
}

#endif // GRAVITY_TREE
