#include "globals.h"
#include "timestep.h"

struct TimeData Time;

void Set_New_Timesteps()
{
	const float a = len3(P[1].Acc);

	float dt = sqrt(2*Param.TimeIntAccuracy*Param.GravSoftening / a);

	Time.Step = dt;

	return ;
}

void Setup_Timesteps()
{
	Time.StepMax = Time.End - Time.Begin;

	return ;
}

bool Time_Is_Up()
{
	if (Time.Current == Time.End)
		return true;

	if (Flag.Endrun)
		return true;

	if (Runtime() >= Param.RuntimeLimit) {

		Flag.WriteRestartFile = true;

		return true;
	}

	return false;
}

bool Time_For_Snapshot()
{
	if (Flag.WriteSnapshot)
		return true;
	
	if (Time.Current == Time.NextSnap) {
	
		Time.NextSnap += Time.BetSnap;

		return true;
	}

	return false;
}

float Timestep(const int ipart)
{
	return Time.Step;
}


