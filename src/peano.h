#define N_PEANO_BITS (sizeof(peanoKey)*CHAR_BIT)
#define N_PEANO_TRIPLETS (N_PEANO_BITS/3)
#define N_SHORT_BITS (sizeof(shortKey)*CHAR_BIT)
#define N_SHORT_TRIPLETS (N_SHORT_BITS/3)

void Sort_Particles_By_Peano_Key();

peanoKey Peano_Key(const double x, const double y, const double z);
peanoKey Reversed_Peano_Key(const double, const double, const double);

shortKey Short_Peano_Key(const float, const float, const float);
shortKey Reversed_Short_Peano_Key(const float, const float, const float);

void test_peanokey();
