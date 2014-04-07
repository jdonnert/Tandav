/* Code units, we have them as macro and constant */

#ifdef UNITS_GADGET_STD

#define LENGTH2CGS 3.0856802e+21
#define MASS2CGS 1.9890000e+43
#define VELOCITY2CGS 1.0e+05
#define TIME2CGS (LENGTH2CGS/VELOCITY2CGS)

#endif // UNITS_GADGET_STD

const struct Units_In_Cgs {
	double Length;
	double Mass;
	double Velocity;
	double Time;
	double Energy;
} Unit;

/* Conversion functions */
double Density_Cgs(const int ipart);
double Temperature_Cgs(const int ipart);
