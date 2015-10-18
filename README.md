tandav - A cosmological simulation code
=======================================

Collisionless N-body dynamics

To Do
-----

1. Main infrastructure
	x macros, globals etc
	x memory management
	x profiling infrastucture
2. I/O
	x Parameter File
	x Reading
	x Writing
	- restart files
3. Cosmology + Units
	x find grav equations
	x transform grav equations to comoving
4. Time integration
	x Cosmo Drift Factor
	x Individual time steps
	x Kepler test problem
x. Peano Hilbert Order 
	x key generation
	x omp sorting
6. Domain Decomposition
	x Find median position of particles in
	  a parallel way to computer CoM w/o
	  PERIODIC
x. Parallel Tree
x. Ewald Summation
	x cube
	x simple
	x tree
	- test
8. Cosmological Sim 
	- test case
	- IDL script for small IC generation
	- growth rate plotting & comparison
9. Multi-grid gravity (Minneapolis & Cray)
	- fast grid definition
10. MPI parallelisation
	- General exchange in domain
	- Engine in tree_accel. Work queues
	  via general stack ?
	- domain decomposition still ok ?
11. Optimisaton
	- single core speed via Cray tools
	- omp speed 
	- mpi & mpi stack stability
	- scalability to 200000 cores
12. Tree on accelerators
	- FMM for run & construction
	- OpenACC general interface
	- communication model
13. MHD via moving mesh
	- here be dragons

