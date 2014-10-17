This is the Style Guide :

* Read: Kernighan & Pike - "The practice of programming"

* IF the code is broader than 79 characters, you are doing it wrong.
	The rationale is that broad code is hard to read. This can usually be
	solved by introducing new variables, which the compiler might later
	optimize away. 
	For example compare :
	
        int myGreatArray = malloc(Task.Npart * N_BINS * 
			sizeof(*my_great_array));
	with 

		size_t nBytes = Task.Npart * N_BINS * sizeof(*myGreatArray);
		int myGreatArray = malloc(nBytes);
	
    Consider a large number of nested loops and conditions. While this
	can not always be avoided, proper modularisation and shorter local variable
	names should make this rarely a problem.

* We format like the Linux kernel - Linus is right, but a tab is just 4 spaces.

* C99

* Variables are initialised when declared. Use the const keyword. Minimise scope.
	
* Global variables have long meaningful names, start with a capital letter.
	Scope should be visible and global variables have to be understandable and 
	unambiguous. You might not want to use these names locally, but define a 
	local variable using the const keyword.  Code-wide variables should be
	embedded in the existing structures, if possible.

* Local variables are short and start with a small letter.

* Modulerise: Every .c file has a corresponding .h header file of the same name. The 
    header file contains the Global functions and variables. 
    These all start with a capital letter and are enclosed in an #ifdef if
    the functionality is switchable. All header file are includes in proto.h
    which itself is contained in globals.h.

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

* Comments are // on the side, /* */ on the line.
	Saves lines, increases readability.

* return ; is not a function.
* The C idiom for infinite loops is for(;;)
* strcpy() is deprecated, always use strncpy(), it's safer.
* All char buffers have size CHARBUFSIZE ! That also helps you to use strncpy.
* Bitmasks are set in hex format : int a = 0x0A;
* mark what you are closing in preprocessor macros : #endif // PMGRID
* inline functions should not contain if () statements depending on function input. 
   Keep them short as well.
* All integers are int for simplicity, if there is no good reason. Otherwise 
  use int64_t, uint32_t etc. long and int are architecture dependent and discouraged.


