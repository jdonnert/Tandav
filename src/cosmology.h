struct CurrentCosmologyInCodeUnits {
	const double Hubble_Constant; // Constants
	const double Omega_Lambda;
	const double Omega_Matter;
	const double Omega_Baryon;
	const double Omega_0;
	const double Omega_Rad;
	const double Rho_Crit0;
	double Hubble_Parameter; // Changing every timestep 
	double Redshift;
	double Expansion_Factor;
	double Critical_Density;
} Cosmo;

double Hubble_Parameter(const double a);
double E_Hubble(const double a);
double Critical_Density(double);

#ifdef COMOVING

void Set_Current_Cosmology();

#else // ! COMOVING

static inline void Set_Current_Cosmology(){};

#endif // ! COMOVING
