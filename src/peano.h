#ifndef PEANO_H
#define PEANO_H

#include "includes.h"
#include "domain.h"

#define N_PEANO_BITS (sizeof(peanoKey)*CHAR_BIT)
#define N_PEANO_TRIPLETS (N_PEANO_BITS/3)
#define N_SHORT_BITS (sizeof(shortKey)*CHAR_BIT)
#define N_SHORT_TRIPLETS (N_SHORT_BITS/3)

#define DELTA_PEANO_BITS (N_PEANO_BITS - N_SHORT_BITS) 

void Sort_Particles_By_Peano_Key();
void Reverse_Peano_Keys();

peanoKey Peano_Key(const Float px, const Float py,const Float pz);
peanoKey Reversed_Peano_Key(const Float px, const Float py, const Float pz);
peanoKey Reverse_Peano_Key(const peanoKey pkey);

shortKey Short_Peano_Key(const Float px, const Float py, const Float pz);
shortKey Reversed_Short_Peano_Key(const Float px, const Float py, 
								  const Float pz);

shortKey Peano2Short_Key(const peanoKey pkey);

void Test_Peanokey();

#endif // PEANO_H
