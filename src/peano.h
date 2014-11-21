#define N_PEANO_TRIPLETS (sizeof(peanoKey)*CHAR_BIT/3)

typedef __uint128_t peanoKey;

void Sort_Particles_By_Peano_Key();
peanoKey Peano_Key(const float, const float, const float);
peanoKey Reversed_Peano_Key(const float, const float, const float);
void test_peanokey();
