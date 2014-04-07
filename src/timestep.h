#ifndef TIMESTEP_H
#define TIMESTEP_H

struct TimeData {
	uint64_t Bin;
	int SnapCounter;
	double StepMax;
	double Begin;
	double Current;
	double End;
	double FirstSnap;
	double BetSnap;
} Time;

void Set_New_Timesteps();
float Timestep(const int);
bool Time_For_Snapshot();
bool Time_Is_Up();

#endif // TIMESTEP_H
