struct Cosmology_Infos { 
	double Omega0;			// Matter content
	double OmegaLambda;		// Dark Energy content
	double OmegaBaryon;		// Baryon content
	double HubbleParam;		// 1% of Hubbles Constant
	double RhoCrit0;		// Critical mass density at z=0
	double Sigma8;			// Norm of matter power-spectrum
	double SpectralIndex;	// Index of matter power-spectrum
	double BaryonFraction;	// Cosmic baryon fraction
} Cosmo; 

double Cosmo_Drift_Factor(double a);
double Cosmo_Kick_Factor(double a);
void Init_Cosmology ();
double Hubble_Function (double a);
double Critical_Density (double a);
