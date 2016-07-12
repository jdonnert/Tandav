#ifndef PERIODIC_H
#define PERIODIC_H

#include "includes.h"

#ifdef PERIODIC

void Periodic_Constrain_Particles_To_Box();
void Periodic_Nearest(Float dr[3]);
void Periodic_Setup();

#else // PERIODIC

static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(Float dr[3]) {};
static inline void Periodic_Setup() {};

#endif //PERIODIC

#endif // PERIODIC_H
// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
