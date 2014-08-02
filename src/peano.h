void test();

typedef uint64_t peanoKey;

void Sort_Particles_By_Peano_Key();
peanoKey Peano_Key(const float x, const float y, const float z, 
		const double *boxsize);
void test_peanokey();
void print_int_bits64(const uint64_t val);
