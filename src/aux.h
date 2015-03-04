#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Warn(...) Warn_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

#define Reallocate_P(...) \
					Reallocate_P_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

void Reallocate_P_Info(const char *, const char *, int, int*, size_t*);

void Assert_Info(const char *, const char *, int, int64_t, const char *, ...);
void Warn_Info(const char *, const char *, int, int64_t, const char *, ...);

/*
 * Helper monkeys
 */

int32_t imin(const int32_t x, const int32_t y); // breaks naming convention :(
int32_t imax(const int32_t x, const int32_t y);
int64_t Imin(const int64_t x, const int64_t y);
int64_t Imax(const int64_t x, const int64_t y);

uint32_t umin(const uint32_t x, const uint32_t y);
uint32_t umax(const uint32_t x, const uint32_t y);
uint64_t Umin(const uint64_t x, const uint64_t y);
uint64_t Umax(const uint64_t x, const uint64_t y);

Float Sign(const Float x);

void Print_Int_Bits(const __uint128_t val, const int length, const int delta);

int Fread(void *restrict data, const size_t size, const size_t nWanted, 
		FILE *stream);
int Fwrite(void *restrict data, const size_t size, const size_t nWrite, 
		FILE *stream);
