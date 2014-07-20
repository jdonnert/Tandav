#ifndef TIMESTEP_H
#define TIMESTEP_H

struct TimeData {
	int SnapCounter;
	float Step; 		// Current Timestep, Min over all Tasks
	double StepMax;
	double StepMin;
	double Begin;
	double Current;		// if COMOVING is set, this is "a"
	double End;
	double FirstSnap;
	double BetSnap;
	double NextSnap;
	int NSnap;
} Time;

void Set_New_Timesteps();
void Setup_Timesteps();
float Timestep(const int);
bool Time_For_Snapshot();
bool Time_Is_Up();

#endif // TIMESTEP_H
