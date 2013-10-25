This is the Style Guide :

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

* We indent like the Linux kernel - Linus is right.
	You can't argue about taste.

* Variables are initialised when declared. 
	Index vars are the obvious exception. This saves lines and the declaration 
	can easily be found using ctags. Furthermore this encourages the use of the 
	const keyword. It also minimises scope, which is always a good thing.
	
* Global variables have long meaningful names, start with a capital letter.
	Obviously the main structure like P, Task and Sim are exceptions, because
	you have to learn those anyway.
	Scope should be visible and global variables have to be understandable and 
	unambiguous. You might not want to use these names locally, but define a 
	local variable using the const keyword. This might even help the compiler.
	In general you should only need to define new global variables on file
	scope ("static X" outside of a function). Code-wide variables should be
	embedded in the existing structures, if possible.

* Local variables are short, hungarian and start with a small letter.
	Again scope should be visible. Short names make it easy to code
	efficiently and stay readable.

* Every .c file has a corresponding .h header file of the same name. The 
    header file contains the Global functions and variables. 
    These all start with a capital letter and are enclosed in an #ifdef if
    the functionality is switchable. All header file are includes in proto.h
    which itself is contained in globals.h.
    The exception to this rule is core of the code, that is cannot be switched 
    off. Here global variables are in globals.h and prototypes are in proto.h.
    This ensures modularisation, whereever possible !

* Constants are macros, have unique long descriptive capitalised names
	There is no elegant alternative to this in pure C. Don't forget to bracket
	everything thats not a one word constant, or division may fail.
	E.g. :
		#define SPEED_OF_LIGHT 3e10
		#define ELECTRON_REST_MASS_ENERGY (me*c*c)

    This will be clunky in your code. Use a constant local variable :
        const double c = SPEED_OF_LIGHT;

    This solves the naming problem and lets you write equations that look like
    math.

* Comments are // on the side, /* */ on the line.
	Saves lines, increases readability.

* return ; is not a function.
* The C idiom for infinite loops is for(;;)
* strcpy() is deprecated, always use strncpy(), it's safer.
* All char buffers have size CHARBUFSIZE ! That also helps you to use strncpy.
* Bitmasks are set in hex format : int a = 0x0A;
* mark what you are closing in preprocessor macros : #endif // PMGRID
* write in C99
* inline functions should not contain if () statements depending on function input. The compiler does not always treat this correctly. Keep them short as well.
* All integers are int for simplicity, if there is no good reason. Otherwise use int64_t, uint64_t etc. long and int are architecture dependent :-(


