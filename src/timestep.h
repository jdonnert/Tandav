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
	uint64_t IntBeg;		// beginning of integer timeline
	uint64_t IntCurrent;	// current point on integer timeline
	uint64_t IntEnd;		// end of integer timeline
	uint64_t IntStep;		// current time step on integer timeline
	uint64_t IntSyncPoint;  // next sync point on integer timeline
	double Step;			// physical time step
	double StepMin;			// smallest physical timestep
	double StepMax;			// largest physical timestep
} Time;

void Set_New_Timesteps();
void Setup_Time_Integration();
double Timebin2Timestep(const int);
double Integer2PhysicalTime(uint64_t);

#endif // TIMESTEP_H
