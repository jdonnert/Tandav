struct TimeData {
	uint64_t Bin;
	int SnapCounter;
	double StepMax;
} Time;

void Set_New_Timesteps();
float Timestep(const int);
bool Time_For_Snapshot();
bool Time_Is_Up();
