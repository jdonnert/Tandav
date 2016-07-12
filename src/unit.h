#ifndef UNIT_H
#define UNIT_H

#include "includes.h"

/* Code unit conversion */

const struct Code_Units {
	double Length;
	double Mass;
	double Velocity;
	double Time;
	double Energy;
	double Density;
} Unit;

void Init_Units();

/* Conversion functions from code to code units */
double Pressure(const int ipart);
double Internal_Energy(const int ipart); // U
double Temperature(const int ipart);

/* Conversion functions to cgs */
double Position_Cgs(const float x);
double Velocity_Cgs(const float v);
double Mass_Cgs(const float mass);
double Density_Cgs(const float rho);
double Number_Density_Cgs(const float rho);
double Pressure_Cgs(const float pres);
double Thermal_Energy_Density_Cgs(const int ipart);

/* Conversion from comoving to physical units */
double Acceleration_Physical(const int ipart);

#endif // UNIT_H

// Copyright (C) 2013 Julius Donnert (donnert@ira.inaf.it)
