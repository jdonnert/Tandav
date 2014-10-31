
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
	double Mean_Density;
} Cosmo;

void Set_Current_Cosmology();
double Redshift();
double Hubble_Parameter(const double a);
double E_Hubble(const double a);
double Critical_Density(double); 
double Mean_Density(double); 
