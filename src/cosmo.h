struct Cosmology_Infos { 
	double Omega0;			// Matter content
	double OmegaLambda;		// Dark Energy content
	double OmegaBaryon;		// Baryon content
	double HubbleParam;		// 1% of Hubbles Constant
	double RhoCrit;			// Critical mass density at z=0
	double Sigma8;			// Norm of matter power-spectrum
	double SpectralIndex;	// Index of matter power-spectrum
	double BaryonFraction;	// Cosmic baryon fraction
} Cosmo; 

void Init_Cosmology ();
float Hubble_Function (float a);
float Critical_Density (float a);
