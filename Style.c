This is the Style Guide :

	“First off, Id suggest printing out a copy of
	 the GNU coding standards, and NOT read it.
	 
	 Burn them, its a great symbolic gesture.”
												~ Linus Torvalds 

* Literature: 
 	- Kernighan & Pike - "The practice of programming"
	- Linus Torvalds on git at google 2007: 
	https://www.youtube.com/watch?v=4XpnKHJAok8
	- http://www.kroah.com/linux/talks/ols_2002_kernel_codingstyle_paper/codingstyle.ps
				

* If code is broader than 79 characters, you are doing it wrong. Except if it
  is a formula. Broad code is hard to read and usually obscures the algorithm.
  This can often be solved by introducing new variables, which the
  compiler might later optimize away, but which explain the algorithm. 
  For example compare :
	
        int my_array = malloc(Task.Npart*N_BINS*sizeof(*my_array));

  with 

		size_t nBytes = Task.Npart * N_BINS * sizeof(*my_array);
		int my_array = malloc(nBytes);
	
  Avoid a large number of nested loops and conditions. Rewrite conditions
  using continue to ease reading, e.g.:

		for (...) {
  
		  	if (A) {
				...
		  	}
		}

		for (...) {
		
			if (!A)
				continue;
			
			...
		}

* Self-explaining code doesn't need many comments. If you modulerize properly
  you will call many static functions whose names will explain most of what
  needs to be known. These function will be optimised away by modern compilers.

* We format like the Linux kernel - Linus is right, but a tab is just 4 spaces.

* C99: Variables are initialised when declared. Use the const keyword for 
  input parameters to avoid bugs. 
  
* Minimize scope! Even declare variables in loop heads like this: 
  for(int i = 0; i < N; i++). This is sometimes even faster in OpenMP.
	
* Global variables should have long meaningful names, start with a capital 
  letter. Scope should be visible and global variables have to be 
  understandable and unambiguous. You might not want to use these names 
  locally, but define a local variable using the const keyword. Code-wide 
  variables should be embedded in the existing structures, if possible. 
  Minimize the use of global variable if reasonably possible.

* Local variables are short and start with a small letter. Don't do this : 
			
			int IamAVeryLongVariableName = 0;
  instead
  			int I_Am_A_Very_Long_Variable_Name = 0;

* If subroutines return multiple values make that clear by declaring them 
  void and return all values by pointer:

  			void my_routine(const int, double *, double *);
  
			my_routine(ipart, &return_var_1, &return_var_2);

* Modulerize: Every .c file has a corresponding .h header file of the same 
  name. The header file contains the Global functions and variables. 
  These all start with a capital letter and are enclosed in an #ifdef if
  the functionality is switchable. All header additional files are includes 
  in proto.h which itself is contained in globals.h.

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
  initialised from the macros.

* return ; is not a function, no brackets ().

* The C idiom for infinite loops is for(;;)

* strcpy() is deprecated, always use strncpy(), it's safer.

* All char buffers have size CHARBUFSIZE ! That also helps you to use strncpy.

* If you fiddle with bits, set constants in hex format : int a = 0x0A;

* Mark what you are closing in preprocessor macros : 

		#ifdef PMGRID
			...
		#endif // PMGRID

* inline functions should not contain if () statements depending on function 
  input. Keep them short as well.

* All integers are int for simplicity, if there is no good reason. Otherwise 
  use int64_t, uint32_t etc. long and int are architecture dependent and
  discouraged.

* Parallelisation is exposed as far up as possible in the call hierachy. 
  Usually that means that all of the MPI, OpenMP and OpenACC calls/directives
  are in the first function of a module, to clearly expose the structure of
  the algorithm. The subroutines then contain the actual physics. Thereby most
  subroutines are NOT thread safe, e.g. all memory function have to be 
  enclosed in #pragma omp single ! This also avoids bugs. 

* All OpenMP globals are public by default and their modification has to be
  protected by single, critical etc. 

* Comments are // on the side, /* */ on the line.
  Saves lines, increases readability.


