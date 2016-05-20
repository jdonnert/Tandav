#ifndef RESTART_FILE_H
#define RESTART_FILE_H

struct Restart_Parameters {
	double Time_Continue;		// hold time if we restart
	double Snap_Counter;
} Restart;

#endif // RESTART_FILE_H
