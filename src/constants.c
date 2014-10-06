#include "globals.h"

#pragma omp threadprivate(Const)
const struct Constants_In_Code_Units Const = {	
	ADIABATIC_INDEX_MONOATOMIC_GAS, // gamma 
	GRAVITATIONAL_CONST/p3(VELOCITY2CGS)/TIME2CGS*MASS2CGS, // Gravity
	HYDROGEN_FRACTION, // xH
	HELIUM_FRACTION // xHe
}; 
