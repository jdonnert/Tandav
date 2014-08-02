#ifndef TIMESTEP_H
#define TIMESTEP_H

typedef uint64_t intime_t; // type of integer time

struct TimeData {
	double Begin;
	double Current;		// if COMOVING is set, this is "a"
	double End;
	int SnapCounter;
	int NSnap;
	double FirstSnap;
	double BetSnap;
	double NextSnap;
	intime_t IntBeg;		// beginning of integer timeline
	intime_t IntCurrent;	// current point on integer timeline
	intime_t IntNext;		// next point on interger timeline
	intime_t IntEnd;		// end of integer timeline
	intime_t IntStep;		// current time step on integer timeline
	intime_t IntFullStep;   // next full step on integer timeline
	double Step;			// physical time step
	double StepMin;			// smallest physical timestep
	double StepMax;			// largest physical timestep
	int MaxActiveBin;		// largest currently active timebin
} Time;

void Set_New_Timesteps();
void Make_Active_Particle_List();

void Setup_Time_Integration();
double Timebin2Timestep(const int);
double Integer2PhysicalTime(intime_t);

#endif // TIMESTEP_H
