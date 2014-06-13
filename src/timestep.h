#ifndef TIMESTEP_H
#define TIMESTEP_H

struct TimeData {
	uint64_t Bin;
	int SnapCounter;
	double StepMax;
	double Begin;
	double Current;
	double Step; 		// Current Timestep
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
