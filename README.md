tandav - A cosmological simulation code
=======================================

Collisionless N-body dynamics

To Do
-----

1. Main infrastructure
	* (DONE) macros, globals etc 
	* (DONE) memory management
	* (DONE) profiling infrastucture
2. I/O
	* (DONE) Parameter File
	* (DONE) Reading
	* (DONE) Writing
	- restart files
3. Cosmology + Units
	* (DONE) find grav equations
	* (DONE) transform grav equations to comoving
4. Time integration
	* (DONE) Cosmo Drift Factor
	* (DONE) Individual time steps
	* (DONE) Kepler test problem
5. Peano Hilbert Order 
	* (DONE) key generation
	* (DONE) omp sorting
6. Domain Decomposition
	* (DONE) Find median position of particles in
	  a parallel way to computer CoM w/o
	  PERIODIC
7. Parallel Tree
 	* (DONE)
8. Ewald Summation
	* (DONE) cube
	* (DONE) simple
	* (DONE) tree
	- test
9. Cosmological Sim 
	* implementation
	- test case
	- IDL script for small IC generation
	- growth rate plotting & comparison
10. Vectorisation 
	- P to struct of arrays
	- Reallocate_P to struct of arrays
	- Adapt code
	- Adapt tree to struct of arrays
11. FMM
	- simple walk
	- fmm_build
	- cache blocking
12. MPI parallelisation
	- General exchange in domain
	- Engine in tree_accel. Work queues
	  via general stack ?
	- domain decomposition still ok ?
13. Optimisaton
	- single core speed via Cray tools
	- omp speed 
	- mpi & mpi stack stability
	- scalability
14. FMM on accelerators
	- OpenACC general interface
	- communication model
15. Multi-grid gravity (Minneapolis & Cray)
	- fast grid definition
16. MHD 
	- derive comoving MHD
	- set thermodynamic model
	- set composition model
	- simple 3D TVD mesh
	- Voronoi generator
		- convex hull

