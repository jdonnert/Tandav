#ifndef ACCEL_H
#define ACCEL_H

#include "includes.h" 
#include "Gravity_Tree/tree.h" 
#include "Gravity/forcetest.h" 
#include "Gravity_FMM/fmm.h"  // <-- add your module .h below

void Compute_Acceleration();
void Safe_Last_Accel();

#ifndef GRAVITY
static inline void Gravity_Acceleration(){};
#endif

#endif // ACCEL_H
