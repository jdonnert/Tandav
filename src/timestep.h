struct TimeData {
	uint64_t Bin;
	int SnapCounter;
	double Current;
} Time;
bool Time_For_Snapshot();
bool Time_Is_Up();
