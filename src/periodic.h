#ifdef PERIODIC

extern void Periodic_Constrain_Particles_To_Box();
Float Periodic_Nearest(const Float dx);
Float Periodic_Nearest_Noncubic(const Float dx, const int i);

#else // !PERIODIC

static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(const Float dx) {};

#endif //PERIODIC

