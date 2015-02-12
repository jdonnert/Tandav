#define N_PEANO_BITS (sizeof(peanoKey)*CHAR_BIT)
#define N_PEANO_TRIPLETS (N_PEANO_BITS/3)
#define N_SHORT_BITS (sizeof(shortKey)*CHAR_BIT)
#define N_SHORT_TRIPLETS (N_SHORT_BITS/3)

void Sort_Particles_By_Peano_Key();

peanoKey Peano_Key(const Float pos[3]);
peanoKey Reversed_Peano_Key(const Float pos[3]);

shortKey Short_Peano_Key(const Float pos[3]);
shortKey Reversed_Short_Peano_Key(const Float pos[3]);

void Test_Peanokey();
