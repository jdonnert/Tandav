#ifndef DRIFT_H
#define DRIFT_H

#include "includes.h"
#include "timestep.h"
#include "periodic.h"
#include "comoving.h"
#include "cosmology.h"
#include "domain.h"
#include "log.h"
#include "Gravity/tree.h"
#include "Gravity/fmm.h"

void Drift_To_Sync_Point();
void Drift_To_Snaptime();
extern double Particle_Drift_Step(const intime_t, const intime_t);

#endif // DRIFT_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
