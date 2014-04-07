void Reallocate_P_Info(const char *, const char *, int, int*, size_t*);

void Assert_Info(const char *, const char *, int, int64_t, const char *, ...);
void Warn_Info(const char *, const char *, int, int64_t, const char *, ...);

float Particle_Mass(const int ipart);
int Particle_Type(const int ipart);
