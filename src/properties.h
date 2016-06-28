#ifndef PROPERTIES_H
#define PROPERTIES_H

#include "includes.h"
#include "particles.h"

void Compute_Current_Simulation_Properties();

struct Simulation_Properties {
	double Total_Mass;			// sum over P.Mass, updated every timestep
	double Center_Of_Mass[3];	// sum of P.Mass*P.Pos
	double Kinetic_Energy;		// sum of 0.5 *P.Mass*P.Vel^2
	double Momentum[3];		    // sum of P.Mass*P.Vel
	double Angular_Momentum[3];	// sum of P.Mass*P.Pos x P.Vel
} Prop;

#endif // PROPERTIES_H
