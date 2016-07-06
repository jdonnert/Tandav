#ifndef GRAVITY_PERIODIC_H
#define GRAVITY_PERIODIC_H

#include "../includes.h"

#if defined(GRAVITY) && defined(PERIODIC)
void Ewald_Correction(const Float dr[3], Float f[3]);
void Gravity_Periodic_Setup();
#else
static inline void Ewald_Correction(const Float dr[3], Float f[3]) {};
static inline void Gravity_Periodic_Setup() {};
#endif // PERIODIC && GRAVITY

#if defined(GRAVITY) && defined(PERIODIC) && defined(GRAVITY_POTENTIAL)
void Ewald_Potential(const Float dr[3], Float p[1]);
#else
static inline void Ewald_Potential(const Float dr[3], Float p[1]) {};
#endif // PERIODIC && GRAVITY && GRAVITY_POTENTIAL

#endif // GRAVITY_PERIODIC_H
