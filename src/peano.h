#define N_PEANO_BITS (sizeof(peanoKey)*CHAR_BIT)
#define N_PEANO_TRIPLETS (N_PEANO_BITS/3)
#define N_SHORT_BITS (sizeof(shortKey)*CHAR_BIT)
#define N_SHORT_TRIPLETS (N_SHORT_BITS/3)

void Sort_Particles_By_Peano_Key();

peanoKey Peano_Key(const Float px, const Float py,const Float pz);
peanoKey Reversed_Peano_Key(const Float px, const Float py,const Float pz);

shortKey Short_Peano_Key(const Float px, const Float py,const Float pz);
shortKey Reversed_Short_Peano_Key(const Float px, const Float py,const Float pz);

void Test_Peanokey();
