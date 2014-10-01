#ifndef TIMESTEP_H

#define TIMESTEP_H

typedef uint64_t intime_t; // type of integer time

struct TimeData {
	double Begin;
	double Current;			// if COMOVING is set, this is "log(a)"
#ifdef COMOVING
	double a;				// expansion factor
#endif
	double End;
	double First_Snap;
	double Bet_Snap;
	double Next_Snap;
	double Step;			// physical time step
	double Step_Min;		// smallest physical timestep
	double Step_Max;		// largest physical timestep
	int Max_Active_Bin;		// largest currently active timebin
	int Step_Counter;
	int Snap_Counter;
	int NSnap;				// expected number of snapshot
} Time;

struct IntergerTimeLine {
	intime_t Beg;		// beginning of integer timeline
	intime_t Current;	// current point on integer timeline
	intime_t Next;		// next point on interger timeline
	intime_t End;		// end of integer timeline
	intime_t Step;		// current time step on integer timeline
	intime_t Full_Step; // next full step on integer timeline
} Int_Time;

void Set_New_Timesteps();

void Setup_Time_Integration();
double Timebin2Timestep(const int);
double Integer2Physical_Time(intime_t);

#endif // TIMESTEP_H
