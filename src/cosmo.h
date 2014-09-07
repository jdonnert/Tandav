struct Cosmology_Infos { 
	double Omega0;			// Matter content
	double Omega_Lambda;		// Dark Energy content
	double Omega_Baryon;		// Baryon content
	double Hubble_Param;		// 1% of Hubbles Constant
	double Rho_Crit0;		// Critical mass density at z=0
	double Sigma8;			// Norm of matter power-spectrum
	double Spectral_Index;	// Index of matter power-spectrum
	double Baryon_Fraction;	// Cosmic baryon fraction
} Cosmo; 

double Cosmo_Drift_Factor(double a);
double Cosmo_Kick_Factor(double a);
void Init_Cosmology ();
double Hubble_Function (double a);
double Critical_Density (double a);
