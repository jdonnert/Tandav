/* This file defines global physical constants in cgsm units
 * Use unique descriptive names */

/* math */
#define PI M_PI
#define DEG2RAD 0.017453293

/* physics */
#define SPEED_OF_LIGHT GSL_CONST_CGSM_SPEED_OF_LIGHT
#define ELECTRON_CHARGE GSL_CONST_CGSM_ELECTRON_CHARGE * GSL_CONST_CGSM_SPEED_OF_LIGHT
#define PLANCK_CONST GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define BOLTZMANN_CONST GSL_CONST_CGSM_BOLTZMANN
#define PROTON_MASS GSL_CONST_CGSM_MASS_PROTON
#define ELECTRON_MASS GSL_CONST_CGSM_MASS_ELECTRON
#define THOMPOSN_CROSS_SECTION GSL_CONST_CGSM_THOMSON_CROSS_SECTION
#define STEPHAN_BOTZMANN_CONST 5.67037321e-8 // [erg/cm^2/s/K^4]
#define FINE_STRUCTURE_CONST GSL_CONST_NUM_FINE_STRUCTURE
#define GRAVITATIONAL_CONST 6.673848e-8 // [cm^3/g/s^2]

/* unit conversions */
#define BARN2CGS GSL_CONST_CGSM_BARN // [cm^2]
#define MEV2CGS 1e6*GSL_CONST_CGSM_ELECTRON_VOLT // [erg]
#define GEV2CGS 1e9*GSL_CONST_CGSM_ELECTRON_VOLT // [erg]
#define MSOL2CGS GSL_CONST_CGSM_SOLAR_MASS // [g]
#define KPC2CGS 3.08567758e21 // [cm]
#define YR2SEC 31556926

/* models */
#define TCMB 2.728 // [K] Temperature of the CMB
#define BCMB = 3.24516e-6 // [G] Magnetic field equivalent of Tcmb
#define H_FRACTION 0.76 // Hydrogen fraction
#define HE_FRACTION 1-H_fraction // Helium fraction
#define UMOL 4.0/(5.0*H_fraction+3.0) // Mean mol. weight in H
#define ADIABATIC_INDEX 5.0/3.0

