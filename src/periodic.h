#ifdef PERIODIC
void Periodic_Constrain_Particles_To_Box();
void Periodic_Nearest(double dr[3]);
void Init_Periodic();

double Boxsize, Boxhalf;
#else
static inline void Periodic_Constrain_Particles_To_Box() {};
static inline void Periodic_Nearest(double dr[3]) {};
static inline void Init_Periodic() {};
#endif //PERIODIC

