#ifndef KICK_H
#define KICK_H

#include "includes.h"
#include "timestep.h"
#include "Gravity/tree.h"
#include "Gravity/fmm.h"

void Kick_First_Halfstep();
void Kick_Second_Halfstep();
extern double Particle_Kick_Step(const intime_t, const intime_t);

#endif // KICK_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
