#ifndef KICK_H
#define KICK_H

#include "includes.h"
#include "timestep.h"
#include "Gravity/tree.h"

void Kick_First_Halfstep();
void Kick_Second_Halfstep();
extern double Particle_Kick_Step(const intime_t, const intime_t);

#endif // KICK_H
