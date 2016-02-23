tandav - A cosmological simulation code
=======================================

Collisionless N-body dynamics

Directories
-----

doc/ 		- Documentation of algorithms
lib/ 		- Plumbing layer to read/write data
src/ 		- All source code
testing/ 	- Library of test problems, including IC generators in IDL and parameter files 

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
	* (DONE) test
9. Cosmological Sim 
	- (DONE) implementation
	- (DONE) test case
	- (DONE) IDL script for IC generation
	- (DONE) Eisenstein & Hu P(k)
	- growth rate plotting & comparison
	- FoF - gadget format ?
	- Halo Output
10. Vectorisation 
	- (DONE) P to struct of arrays
	- (DONE) Reallocate_P to struct of arrays
	- (DONE) Adapt code
	- (DONE) Adapt tree to struct of arrays
	- (DONE) Tree accel vectorize particle loops
	- Particle loop to vector blocks
11. FMM
	- simple walk, same Tree, measure speed
	- fmm_build & cache blocking
	- vectorization
12. MPI parallelisation
	- General exchange in domain
	- Engine in accel. Work queues via general stack ?
	- domain decomposition still ok ?
13. Optimisaton
	- single core speed via Cray tools
	- omp speed 
	- mpi & mpi stack stability
	- scalability
16. MHD 
	- derive comoving MHD
	- set thermodynamic model
	- set composition model
	- simple 3D TVD mesh (->Tom)
	- Voronoi generator
		- convex hull - MPI

