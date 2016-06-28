#ifndef PROFILE_H
#define PROFILE_H

#include  "includes.h"

void Init_Profiler();
void Finish_Profiler();
void Profile_Info(const char* file, const char* func, const int line, 
		const char *name);
void Profile_Report(FILE *);
void Profile_Report_Last(FILE *);
void Write_Logs();
double Runtime();

#endif // PROFILE_H
