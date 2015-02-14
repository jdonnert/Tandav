void Drift_To_Sync_Point();
void Drift_To_Snaptime();

double Particle_Drift_Step(const int ipart, const double time_next);

#ifdef PERIODIC
void Periodic_Constrain_Particles_To_Box();
#else // !PERIODIC
static inline void Periodic_Constrain_Particles_To_Box() {};
#endif //PERIODIC
