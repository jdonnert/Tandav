This is the Style Guide for Tandav:

	“First off, I'd suggest printing out a copy of
	 the GNU coding standards, and NOT read it.

	 Burn them, its a great symbolic gesture.”
												~ Linus Torvalds

* Literature:
	- Kernighan & Pike - "The practice of programming"
	- http://youtu.be/LLBrBBImJt4
	- Linus Torvalds on git 2007:
	  https://www.youtube.com/watch?v=4XpnKHJAok8
	- http://www.kroah.com/linux/talks/ols_2002_kernel_codingstyle_paper/codingstyle.ps


* We format like the Linux kernel, but a tab is just 4 spaces 
  (helps with long formulas).

* If code is broader than 80 characters, you are doing it wrong. This is basic 
  typography: reading text becomes harder, when it is broader than 90 characters. 
  Broad code is not only hard to read, but usually also obscures the algorithm.
  This can often be solved by introducing new variables, which the
  compiler might later optimize away, but which explain the algorithm.
  For example compare :

        int *my_array = malloc(3*Task.Npart[1]*N_BINS*sizeof(*my_array));

  with

		size_t nBytes = 3 * Task.Npart[3] * N_BINS * sizeof(int);
		int *my_array = malloc(nBytes);

* Self-explaining code doesn't need many comments, you should use functions. 
  If you modulerize properly you will call many static functions whose names 
  will explain most of what needs to be known. These function will be
  optimised away by modern compilers. Across files, -flto or -ipo switches do
  the same.

* Write short functions, whose name is a description of what you are doing.
  No comment necessary to explain how something is happening, only what.
  
* Avoid a large number of nested loops and conditions. Rewrite conditions
  using continue to ease reading and reduce indentation level, e.g.:

		for (...) {

			if (A) {
				...	lots of code
			} // if A
		}

		for (...) {

			if (!A)
				continue;

			lots of code ...

		}

  An exception for this rule are hot loops, because this approach does not
  vectorize !

* The naming scheme of the modules should be consistent on the Makefile, file 
  and function level:

			#DEFINE GRAVITY_TREE
			src/Gravity/tree_build.c
			Gravity_Tree_Build();

* C99: Variables are initialised when declared. Use the const keyword for
  input parameters of functions, to make clear what is returned.

* In general its a good idea to avoid the pointer picture where possible.
  E.g. if you declare pointers as function parameters and you know their size
  beforehand, tell the reader and the compiler: 

			static void find_domain_center(double Center_Out[3]);

  instead of

			static void find_domain_center(double *Center_Out);

  Now the compiler can in principle check for out of bounds accesses.

* Minimize scope! Even declare variables in loop heads like this:
  for(int i = 0; i < N; i++). This is sometimes even faster in OpenMP 
  and helps the parallelisation in general.

* Global variables should have long meaningful names, start with a capital
  letters. Scope should be visible and global variables have to be
  understandable and unambiguous. This also helps with OpenMP race conditions,
  global variables should be touched only inside a #pragma omp single region.
  You might not want to use long names locally, but define a local variable
  using the const keyword.

  Globally:
            const double Speed_Of_Light = SPEED_OF_LIGHT / Unit.Velocity;
  Locally: 
            const double c = Speed_Of_Light;
  
  Code-wide variables should be embedded in the existing structures, if 
  possible. Minimize the use of global variables if reasonably possible.

* Local variables are short and start with a small letter. Don't do this :

			int IamAVeryLongVariableName = 0;
  
  instead use underscores to increase readability

			int I_Am_A_Very_Long_Variable_Name = 0;

  See Greg KH talk on the kernel

* If subroutines return multiple values make that clear by declaring them
  void and return all values by pointer:

			void my_routine(const int ipart, double * result1, double * result2);

			my_routine(ipart, &return_var_1, &return_var_2);

* Modulerize: Every .c file/module has a corresponding .h header file of the same 
  name. The header file contains the Global functions and variables. 
  These all start with a capital letter and are enclosed in an #ifdef if
  the functionality is switchable. All additional header files are includes 
  in proto.h , which itself is contained in globals.h.

* Constants in CGSM are macros, have unique long descriptive capitalised names
  There is no elegant alternative to this in pure C. Don't forget to bracket
  everything that's not a one word constant, or division may fail.
  E.g. :

		#define SPEED_OF_LIGHT 3e10
		#define ELECTRON_REST_MASS_ENERGY (me*c*c)

  Then use a constant local variable :

        const double c = SPEED_OF_LIGHT;

  This solves the naming problem and lets you write equations that look like
  math.
  Constants in Code units can be found in the Const structure, which is 
  initialised from the macros at compile time. This way its also clear, which
  units a constant carries.

* return ; is not a function, no brackets ().

* goto is OK only (!) if you skip forward inside one function. Then it is 
  encouraged, because it actually simplifies the code. 

* The C idiom for infinite loops is for(;;).

* strcpy() is deprecated, always use strncpy(), it's safer.

* All char buffers have size CHARBUFSIZE ! That also helps you to use strncpy.

* If you fiddle with bits, set constants in hex format : int mask = 0x0A;

* Mark what you are closing in preprocessor macros :

		#ifdef PMGRID
			...
		#endif // PMGRID

  You never know what someone else is going to squeeze into your define later
  so the endif might appear pages down. A negative example is Gadget-3's hydra.c

* Minimize #ifdefs in C code. Write a function and an empty 
  "inline void F(){};" prototype in the header file. Start the function 
  name with the macro name.  See the handling of 
  Periodic_Constrain_Particles_To_Box() in drift.[ch]

* Avoid stacking #ifdef, it becomes unreadable too quickly (hydra.c anybody?) 
  Check Gravity/gravity.h to see how to do it.

* Default integers should be simply int. If you need more bits, unsigned etc, 
  exclusively use int64_t, uint32_t etc. Standard long and unsigned int are 
  architecture dependent. Array & Malloc sizes should be size_t, pointer 
  positions ptrdiff_t.

* Parallelisation is exposed as far up as possible in the call hierachy.
  Usually that means that most of the MPI & OpenMP calls/directives
  are in the first function of a module, to clearly expose the structure of
  the algorithm. The subroutines then contain the actual physics. Thereby most
  subroutines are NOT thread safe, e.g. often memory functions have to be
  enclosed in #pragma omp single ! This avoids bugs as a #pragma critical 
  inside a #pragma single will lead to a segfault and erratic behaviour.
  In some cases, MPI communication is still done in its own function.

* All OpenMP globals are public by default and their modification has to be
  protected by single, critical, etc to avoid race conditions. Sig is private.

* Comments are // on the side, /* */ on the line. 
  Saves lines, increases readability.

* Modules have a common structure: A module X can have Init_X(), Setup_X()
  and Finish_X() functions, defined in the module file and called in init.c 
  setup.c and finish.c. This way you can execute memory allocation etc
  at various stages in code, in particular outside of the omp parallel region.
  This is a C way of writing object oriented code.

* We use the structures of arrays approach, not the array of structures 
  approach. The simple reason is aligned memory access for vectorization. 
  I.e. a loop like this
  			
  			for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
				P.Pos[0][ipart] += P.Vel[0][ipart] * dt;

	vectorizes. This one doesn't :

  			for (int ipart = 0; ipart < Task.Npart_Total; ipart++)
				P[ipart].Pos[0] += P[ipart].Vel[0] * dt;

	The difference in speed is easily a factor of ten !

* We (will) use thread multiple MPI and RMA to overlap work and communication.
