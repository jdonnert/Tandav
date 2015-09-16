/* 
 * This file defines global physical constants in cgsm units
 * Use unique descriptive names, all Capital letters 
 */

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>

#define PI M_PI
#define SPEED_OF_LIGHT GSL_CONST_CGSM_SPEED_OF_LIGHT
#define PLANCK_CONST GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define BOLTZMANN_CONST GSL_CONST_CGSM_BOLTZMANN
#define ELECTRON_CHARGE (GSL_CONST_CGSM_ELECTRON_CHARGE * SPEED_OF_LIGHT)
#define ELECTRON_MASS GSL_CONST_CGSM_MASS_ELECTRON
#define PROTON_MASS GSL_CONST_CGSM_MASS_PROTON
#define THOMPSON_CROSS_SECTION GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define FINE_STRUCTURE_CONST GSL_CONST_NUM_FINE_STRUCTURE
#define STEPHAN_BOLTZMANN_CONST 5.67037321e-8 // [erg/cm^2/s/K^4]
#define GRAVITATIONAL_CONST 6.673848e-8 // [cm^3/g/s^2]

/* 
 * Unit conversions to cgs
 */

#define BARN2CGS GSL_CONST_CGSM_BARN // [cm^2]
#define MEV2CGS (1e6*GSL_CONST_CGSM_ELECTRON_VOLT) // [erg]
#define GEV2CGS (1e9*GSL_CONST_CGSM_ELECTRON_VOLT) // [erg]
#define MSOL2CGS GSL_CONST_CGSM_SOLAR_MASS // [g]
#define KPC2CGS 3.08567758e21 // [cm]
#define MPC2CGS (1e3*KPC2CGS) // [cm]
#define KM2CGS 1e5 // [cm]
#define YR2CGS 31556926.0 // [sec], seconds in a year

/* 
 * Model constants 
 */

#define TCMB 2.728		// [K] Temperature of the CMB
#define BCMB 3.24516e-6	// [G] Magnetic field equivalent of Tcmb
#define HYDROGEN_FRACTION 0.76		// Std hydrogen mass fraction
#define HELIUM_FRACTION (1-HYDROGEN_FRACTION)	// Helium mass fraction
#define UMOL (4.0/(5.0*HYDROGEN_FRACTION+3.0)) // Mean mol. weight in H
#define ADIABATIC_INDEX_MONOATOMIC_GAS (5.0/3.0)

const double Pi, Deg2Rad, Sqrt3, Sqrt2; // math

struct Constants_In_Code_Units {
	const double Adiabatic_Index;
	const double Gravity;
	const double Hydrogen_Fraction;
	const double Helium_Fraction;
} Const ;
