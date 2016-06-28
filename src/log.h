#ifndef LOG_H
#define LOG_H

#include "includes.h"
#include "timestep.h"
#include "properties.h"
#include "IO/parameter_file.h"

/* 
 * This provides a logging infrastructure
 */

void Write_Logs();
void Init_Logs();
void Finish_Logs();

#endif // LOG_H
