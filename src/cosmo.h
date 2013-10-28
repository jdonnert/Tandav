struct Cosmology_Infos { 
	double Omega0;				// Matter content
	double OmegaLambda;			// Dark Energy content
	double OmegaBaryon;			// Baryon content
	double HubbleParam;			// 1% of Hubbles Constant
	double RhoCrit;				// Critical mass density
	double Sigma8;				// Norm of matter power-spectrum
	double SpectralIndex;		// Index of matter power-spectrum
	double BaryonFraction;		// Cosmic baryon fraction
} Cosmo; 

void Fill_Comoving_Factor_Tables();
