#ifdef COMOVING
void Setup_Comoving();
void Finish_Comoving();
double Particle_Drift_Step(const int ipart, const double time_next);
double Particle_Kick_Step(const int ipart, const double time_next);
double Comoving_VelDisp_Timestep_Constraint();
#else
inline void Setup_Comoving() {};
inline void Finish_Comoving() {};
inline double Comoving_VelDisp_Timestep_Constraint(double dt) {return dt;};
#endif // COMOVING
