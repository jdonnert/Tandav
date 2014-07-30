#define Assert(...) Assert_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define Warn(...) Warn_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

#define Reallocate_P(...) Reallocate_P_Info(__func__, __FILE__, __LINE__, __VA_ARGS__)

void Reallocate_P_Info(const char *, const char *, int, int*, size_t*);

void Assert_Info(const char *, const char *, int, int64_t, const char *, ...);
void Warn_Info(const char *, const char *, int, int64_t, const char *, ...);

float Particle_Mass(const int ipart);
int Particle_Type(const int ipart);

int64_t Imin(const int64_t x, const int64_t y);
int64_t Imax(const int64_t x, const int64_t y);
uint64_t Umin(const uint64_t x, const uint64_t y);
uint64_t Umax(const uint64_t x, const uint64_t y);
