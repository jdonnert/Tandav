void Drift_To_Sync_Point();
void Drift_To_Snaptime();

#ifdef PERIODIC
void Constrain_Particles_To_Box();
#else // !PERIODIC
static inline void Constrain_Particles_To_Box() {}
#endif //PERIODIC
