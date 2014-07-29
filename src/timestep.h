#ifndef TIMESTEP_H
#define TIMESTEP_H

struct TimeData {
	double Begin;
	double Current;		// if COMOVING is set, this is "a"
	double End;
	int SnapCounter;
	int NSnap;
	double FirstSnap;
	double BetSnap;
	double NextSnap;
	uint64_t IntBeg;		// == 0	
	uint64_t IntCurrent;	// current point on integer timeline
	uint64_t IntEnd;		// == (1ULL << 63)
	uint64_t IntStep;		// current time step on integer timeline
	uint64_t IntSyncPoint;  // next sync point on integer timeline
	double Step;			// physical time step
	double StepMin;			// smallest physical timestep
	double StepMax;			// largest physical timestep
} Time;

void Set_New_Timesteps();
void Setup_Timesteps();
double Timebin2Timestep(const int);
double Integer2PhysicalTime(uint64_t);
bool Time_For_Snapshot();
bool Time_Is_Up();

#endif // TIMESTEP_H
