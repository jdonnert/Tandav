#include "../globals.h"
#include "gravity.h"
#include "../domain.h"

#ifdef GRAVITY_TREE


static struct GravityDataForExport {
	int Target_Task;
	Float Pos[3];
	Float Acc;
} DataOut;

static struct GravityDataForImport {
	int Target_Part;
	Float Acc[3];
	Float Cost;
} DataIn;

static void gravity_tree_walk_local(const int ipart, Float*, Float*);

static void interact(const Float,const Float *,const Float,Float *,Float *);

static inline Float node_size(const int node);

static inline int level(const int node)
{
	const uint32_t bitmask = 0x3F;

	return Tree[node].Bitfield & bitmask; // return but 0-5
}

void Gravity_Tree_Acceleration()
{
	//Accel_Gravity_Simple();

	Profile("Grav Tree Walk");

	double max_rel_err = 0;
	double mean_err = 0;
	int worst_part = -1;
	int cnt = 0;

	#pragma omp for
	for (int i = 0; i < NActive_Particles; i++) {
			
		int ipart = Active_Particle_List[i];

		Float grav_accel[3] = { 0 }, pot = 0, 

		gravity_tree_walk_local(ipart, grav_accel, &pot);

		/*double rel_err = fabs(ALENGTH3(grav_accel)-ALENGTH3(P[ipart].Acc))
			/ ALENGTH3(P[ipart].Acc);

		if (rel_err > max_rel_err) {
		
			worst_part = ipart;
			max_rel_err = rel_err;
		}

		mean_err += rel_err;

		if (rel_err > 0.01)
			cnt++;

		//printf("%d %d %g %g %g %g\n", ipart, P[ipart].ID, (ALENGTH3(grav_accel)-ALENGTH3(P[ipart].Acc)) / ALENGTH3(P[ipart].Acc), P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2]);

if (rel_err > 0.2)
printf("ipart %d, rel err %g | %g %g %g | %g %g %g| %g %g %g |%g %g %g \n",
				ipart, rel_err,
				grav_accel[0],grav_accel[1],grav_accel[2],
				P[ipart].Acc[0],P[ipart].Acc[1],P[ipart].Acc[2],
				(grav_accel[0] - P[ipart].Acc[0])/grav_accel[0],
				(grav_accel[1] - P[ipart].Acc[1])/grav_accel[1],
				(grav_accel[2] - P[ipart].Acc[2])/grav_accel[2], 
				 P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2]); */

		P[ipart].Acc[0]  += grav_accel[0];
		P[ipart].Acc[1]  += grav_accel[1];
		P[ipart].Acc[2]  += grav_accel[2];

#ifdef OUTPUT_PARTIAL_ACCELERATIONS
		P[ipart].Grav_Acc[0] = grav_accel[0];
		P[ipart].Grav_Acc[1] = grav_accel[1];
		P[ipart].Grav_Acc[2] = grav_accel[2];
#endif 

#ifdef GRAVITY_POTENTIAL
		P[ipart].Grav_Pot = pot;
#endif

	} // ipart
	
//printf("Largest Err %g at %d, %d above threshold, mean err %g \n", 	max_rel_err, worst_part, cnt, mean_err / Task.Npart_Total);

	Profile("Grav Tree Walk");

	return ;
}
	

/*
 * This function walks the local tree and computes the gravitational 
 * acceleration.
 */

static void gravity_tree_walk_local(const int ipart, Float* accel)
{
	int node = 1;

	const Float fac = ALENGTH3(P[ipart].Acc)/Const.Gravity*TREE_OPEN_PARAM_REL;

	while (Tree[node].DNext != 0) {
		
		if (Tree[node].DNext < 0) { // encountered particle bundle

			int first = -(Tree[node].DNext + 1); // part index is offset by 1
			int last = first + Tree[node].Npart;

			for (int jpart = first; jpart < last; jpart++ ) {

				if (jpart == ipart)
					continue;
			
				Float dr[3] = { P[ipart].Pos[0] - P[jpart].Pos[0],	
					   		    P[ipart].Pos[1] - P[jpart].Pos[1], 
					            P[ipart].Pos[2] - P[jpart].Pos[2] };

				Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

				Float mpart = P[jpart].Mass;

				interact(mpart, dr, r2, accel); 
			}

			node++;
			
			continue;
		}
	
		Float dr[3] = { P[ipart].Pos[0] - Tree[node].CoM[0],	
					    P[ipart].Pos[1] - Tree[node].CoM[1], 
					    P[ipart].Pos[2] - Tree[node].CoM[2] };

		Float r2 = p2(dr[0]) + p2(dr[1]) + p2(dr[2]);

		Float nMass = Tree[node].Mass;

		Float nSize = node_size(node); // now check opening criteria

		if (fac == 0) {

			if (nSize*nSize > r2 * TREE_OPEN_PARAM_BH) { // BH criterion

				node++; // open

				continue;
			}

		} else {
		
			if (nMass*nSize*nSize > r2*r2 * fac) { // relative criterion

				node++; 

				continue;
			}
		
		}

		//Float dx = fabs(P[ipart].Pos[0] - Tree[node].Pos[0]);

		if ( dr[0] < 0.6 * nSize) {

			//Float dy = fabs(P[ipart].Pos[1] - Tree[node].Pos[1]);
			
			if ( dr[1] < 0.6 * nSize) {
			
				//Float dz = fabs(P[ipart].Pos[2] - Tree[node].Pos[2]);
				
				if ( dr[2] < 0.6 * nSize) {

					node++; // part inside node

					continue;
				}
			}
		}
		
		interact(nMass, dr, r2, accel); // use node
		
		node += fmax(1, Tree[node].DNext);
		
	} // while

	return ;
}

static const Float h_grav = GRAV_SOFTENING / 3.0; // Plummer equiv softening

static void interact(const Float mass, const Float dr[3], const Float r2, 
		Float accel[3], Float pot)
{
	const Float r = sqrt(r2);
	
	Float rinv = 1/r;

#ifdef GRAVITY_POTENTIAL
	Float rinv_pot = rinv;
#endif

	if (r < h_grav) {
	
		Float u = r/h_grav;
		Float u2 = u*u;
		Float u3 = u2*u;
			
		rinv = sqrt(14*u- 84*u3 + 140*u2*u2 - 90*u2*u3 + 21*u3*u3) / h_grav;

#ifdef GRAVITY_POTENTIAL
		rinv_pot = (7*u2-21*u2*u2+28*u3*u2-15*u3*u3+u3*u3*u*8-3)/h;
#endif
	} 
	
	Float acc_mag = Const.Gravity * mass * p2(rinv);

	accel[0] += -acc_mag * dr[0] * rinv;
	accel[1] += -acc_mag * dr[1] * rinv;
	accel[2] += -acc_mag * dr[2] * rinv;

#ifdef GRAVITY_POTENTIAL
	P[ipart].Grav_Pot += -Const.Gravity * P[jpart].Mass *rinv_pot;
#endif

	return ;
}

static inline Float node_size(const int node)
{
	int lvl = Tree[node].Bitfield & 0x3F; // level

	return Domain.Size / ((Float) (1UL << lvl)); // Domain.Size/2^level
}

#endif // GRAVITY_TREE
