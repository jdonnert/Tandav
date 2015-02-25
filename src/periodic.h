#ifdef PERIODIC

extern void Periodic_Constrain_Particles_To_Box();
extern inline Float Periodic_Nearest(const Float dx);

#else // !PERIODIC

static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(const Float dx) {};

#endif //PERIODIC

