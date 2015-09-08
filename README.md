tandav
======

A cosmological simulation code

Constants
---------

Constants usually come in up to three units: Code units, CGSM and non-standard. 
	
A good example is the Hubble constant usually expressed as H0 = 100 km/Mpc/s,
in CGSM: 3.2407765e-18 1/s, in standard code units: 0.1

Inside the code all constants/conversion factors IN CODE UNITS are stored in 
structures (Unit, Cosmo, Const ...), which are initialised at compile time or during start-up. We do this so modules can make constants time dependent, if they are switched on, e.g. the hydrogen fraction. In this case, an update function can be called once per timestep that set the new value of the constant. Similar to `Set_Current_Cosmology()`.

All other constants are defined as macros, hence represent immutable standard values. This includes the constants defining the code units, and the cosmological model. If you want to change these, you have to introduce them into the structures.

This way it is also obvious that all constants in the structures are in code units. The macro constants can be in CGSM or non-standard (-> H0), here you have to be careful.
