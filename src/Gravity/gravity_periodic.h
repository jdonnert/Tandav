#if defined(GRAVITY) && defined(PERIODIC)
void Ewald_Correction(const double dr[3], Float f[3]);
void Gravity_Periodic_Init();
double Boxsize, Boxhalf; // now the box _has_ to be cubic
#else
static inline void Ewald_Correction(const double dr[3], Float f[3]) {};
static inline void Gravity_Periodic_Init() {};
#endif // PERIODIC && GRAVITY

#if defined(GRAVITY) && defined(PERIODIC) && defined(GRAVITY_POTENTIAL)
void Ewald_Potential(const double dr[3], Float p[1]);
#else
static inline void Ewald_Potential(const double dr[3], Float p[1]) {};
#endif // PERIODIC && GRAVITY && GRAVITY_POTENTIAL
