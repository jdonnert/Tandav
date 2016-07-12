#ifndef FORCETEST_H
#define FORCETEST_H

#include "../includes.h"

#if defined(GRAVITY) && defined(GRAVITY_FORCETEST)
void Gravity_Forcetest();
#else
static inline void Gravity_Forcetest() {};
#endif // GRAVITY_FORCETEST

#endif // GRAVITY_SIMPLE_H





// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
