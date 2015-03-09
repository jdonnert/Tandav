#ifdef PERIODIC

extern void Periodic_Constrain_Particles_To_Box();
void Periodic_Nearest(Float dr[3]);
Float Periodic_Nearest_Noncubic(const Float dx, const int i);

#else // !PERIODIC

static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(Float dr[3]) {};
static inline void Periodic_Nearest_Noncubic(const Float dx, const int i) {};

#endif //PERIODIC

