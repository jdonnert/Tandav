#ifndef ACCEL_H
#define ACCEL_H

#include "includes.h" 
#include "particles.h" 
#include "Gravity/tree.h" 
#include "Gravity/forcetest.h" 
#include "Gravity/fmm.h"  // <-- add your module .h below

void Compute_Acceleration();
void Safe_Last_Accel();

#ifndef GRAVITY
static inline void Gravity_Acceleration(){};
#endif

#endif // ACCEL_H
Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
