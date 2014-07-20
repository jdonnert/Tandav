#include "globals.h"
#include "kick.h"
#include "force.h"
#include "timestep.h"

/* update velocities from forces */
void Kick_Halfstep() 
{
#ifdef COMOVING
	const float driftfac = Cosmo_Kick_Factor(Sim.CurrentTime);
#else
	const float driftfac = 1;
#endif // COMOVING
	
#pragma omp parallel for 
	for (int ipart = 0; ipart < Task.NpartTotal; ipart++) {

		float mpart = P[ipart].Mass;

		float force[3] = { 0 };
		
		float dt = Time.Step;//Timestep(ipart);

		P[ipart].Vel[0] += 0.5 * dt * P[ipart].Force[0]/mpart * driftfac;
		P[ipart].Vel[1] += 0.5 * dt * P[ipart].Force[1]/mpart * driftfac; 
		P[ipart].Vel[2] += 0.5 * dt * P[ipart].Force[2]/mpart * driftfac;
	}

	return ;
}

