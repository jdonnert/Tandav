#ifndef TIMESTEP_H

#define TIMESTEP_H

struct TimeData {
	double Begin;			// if COMOVING is set, these are expansion fact a
	double Current;
	double Next;
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

struct IntegerTimeLine {
	intime_t Beg;			// beginning of integer timeline
	intime_t Current;		// current point on integer timeline
	intime_t Next;			// next point on interger timeline
	intime_t End;			// end of integer timeline
	intime_t Step;			// current time step on integer timeline
	intime_t Next_Sync_Point;// next full step on integer timeline
} Int_Time;

void Set_New_Timesteps();
void Make_Active_Particle_List();
void Make_Active_Particle_Vectors();

void Setup_Time_Integration();
double Timebin2Timestep(const int);
double Integer2Physical_Time(const intime_t);
double Integer_Time2Integration_Time(const intime_t);
intime_t Integration_Time2Integer_Time(const double);

#endif // TIMESTEP_H
