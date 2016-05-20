#ifndef GRAVITY_SIMPLE_H
#define GRAVITY_SIMPLE_H

#include "../includes.h"

#if defined(GRAVITY) && defined(GRAVITY_FORCETEST)
void Gravity_Simple_Accel();
#else
static inline void Gravity_Simple_Accel() {};
#endif // GRAVITY_FORCETEST

#endif // GRAVITY_SIMPLE_H



