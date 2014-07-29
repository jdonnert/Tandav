#include "globals.h"
#include "timestep.h"

static void setup_timeline();

void Setup() 
{
	setup_timeline();

	Setup_Timesteps();

#ifdef COMOVING
	Setup_Cosmology();

	Setup_Comoving();
#endif // COMOVING

	return ;
}

static void setup_timeline()
{
	Time.NextSnap = Time.FirstSnap;

	Time.NSnap = (Time.End - Time.Begin) / Time.BetSnap;

	Assert(Time.NSnap > 0,"Timeline does not seem to produce any outputs");

	rprintf("Simulation timeline: start time = %g, end time = %g,"
			" Number of Snaps = %d \n\n", 
			Time.Begin, Time.End, Time.NSnap);

	return ;
}
