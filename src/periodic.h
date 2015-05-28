#ifdef PERIODIC
void Periodic_Constrain_Particles_To_Box();
void Periodic_Nearest(Float dr[3]);
#else
static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(Float dr[3]) {};
#endif //PERIODIC

