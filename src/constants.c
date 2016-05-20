#include "constants.h"

struct Constants_In_Code_Units Const = {
	ADIABATIC_INDEX_MONOATOMIC_GAS,
	GRAVITATIONAL_CONST/p3(VELOCITY2CGS)/(LENGTH2CGS/VELOCITY2CGS)*MASS2CGS,
	HYDROGEN_FRACTION,
	HELIUM_FRACTION
};

const double Pi = PI;
const double Deg2Rad = PI / 180.0;
const double Sqrt2 = 1.4142135623730951454746218587;
const double Sqrt3 = 1.7320508075688771931766041234;

void Init_Constants()
{
	printf("Constants (internal units):\n"
			"  Gravity = %g \n"
			"  Adiabatic_Index = %g \n"
			"  Hydrogen_Fraction = %g \n"
			"  Helium_Fraction = %g \n\n", 
			Const.Gravity, Const.Adiabatic_Index, Const.Hydrogen_Fraction,
			Const.Helium_Fraction);

	return ;
}
