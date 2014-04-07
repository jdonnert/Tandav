#include "globals.h"
#include "timestep.h"

void Set_New_Timesteps()
{
	return ;
}

void Init_Timesteps()
{
	Time.StepMax = 1;

	return ;
}

bool Time_Is_Up()
{

	return 0;
}

bool Time_For_Snapshot()
{

	return 0;
}

float Timestep(const int ipart)
{
	//return 62 - imax(0 , ceil(log2(Time.StepMax/dt))) - 2;
	return 1;
}


