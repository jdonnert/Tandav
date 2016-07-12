#ifndef COMOV_H
#define COMOV_H

#include "includes.h"
#include "cosmology.h"

#ifdef COMOVING
extern void Comoving_Setup();
extern void Finish_Comoving();
extern double Comoving_VelDisp_Timestep_Constraint();
#else
inline void Comoving_Setup() {};
inline void Finish_Comoving() {};
inline double Comoving_VelDisp_Timestep_Constraint(double dt) {return dt;};
#endif // COMOVING

#endif // COMOV_H
Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
