#ifdef COMOVING
void Setup_Comoving();
void Finish_Comoving();
double Particle_Drift_Step(const int ipart, const double time_next);
double Particle_Kick_Step(const int ipart);
#else
inline void Setup_Comoving() {};
inline void Finish_Comoving() {};
#endif // COMOVING
