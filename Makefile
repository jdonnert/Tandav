# Compilation parameters in 'Config' file

SHELL = /bin/bash # This should always be present in a Makefile

ifndef SYSTYPE # set this in your ~/.bashrc or ~/.tcshrc
SYSTYPE := ${shell hostname}
endif

CC	 	 = mpicc
OPTIMIZE = -Wall -g -O2 
MPI_INCL = $(CPPFLAGS)
MPI_LIBS = $(LDFLAGS)
GSL_INCL =
GSL_LIBS = 
FFT_INCL =
FFT_LIBS =

# machine specifics
ifeq ($(SYSTYPE),DARWIN)
CC       =  mpicc
OPTIMIZE = -O3 -g -Wall -lmpich -mtune=native -march=corei7 #-ftree-vectorizer-verbose=2
MPI_LIBS = 
MPI_INCL = 
GSL_INCL =  
GSL_LIBS = 
FFT_INCL =
FFT_LIBS =
endif

ifeq ($(SYSTYPE),mach64.ira.inaf.it)
CC       =  mpicc
OPTIMIZE =  -g -O2  -march=bdver1 -mtune=native -mprefer-avx128 -mieee-fp  \
			-minline-all-stringops -fprefetch-loop-arrays --param prefetch-latency=300 -funroll-all-loops 
MPI_LIBS = -L/homes/donnert/Libs/lib
MPI_INCL = -I/homes/donnert/Libs/include
GSL_INCL =  
GSL_LIBS = 
FFT_INCL =
FFT_LIBS =
endif

ifeq ($(SYSTYPE),getorin.ira.inaf.it)
CC       =  mpicc
OPTIMIZE =  -Wall -g -O3 -openmp  -finline -finline-functions \
			-funroll-loops  -xhost  -mkl  -ipo 
MPI_LIBS = -lmpich -L/homes/donnert/Libs/lib 
MPI_INCL = -I/homes/donnert/Libs/include 
GSL_INCL =  
GSL_LIBS = -L/opt/intel/composer_xe_2013_sp1.2.144/lib/
FFT_INCL =
FFT_LIBS =
endif

# end systypes

EXEC = Tandav

SRCDIR = src/

SRCFILES = main.c aux.c cosmology.c domain.c update.c print_settings.c drift.c \
		init.c kick.c setup.c timestep.c unit.c memory.c profile.c \
		sort.c finish.c peano.c accel.c constants.c log.c signal.c comov.c \
	   	IO/io.c IO/read_snapshot.c IO/write_snapshot.c IO/rw_parameter_file.c \
		IO/write_restart_file.c IO/read_restart_file.c  \
		Gravity/gravity_simple.c Gravity/gravity_tree.c \

INCLFILES = config.h globals.h cosmology.h unit.h aux.h macro.h proto.h \
	    memory.h profile.h IO/io.h constants.h kick.h setup.h update.h \
		drift.h timestep.h peano.h accel.h log.h signal.h  comov.h \
		Gravity/gravity.h \
		../Makefile ../Config

OBJFILES = $(SRCFILES:.c=.o)

OBJS = $(addprefix $(SRCDIR),$(OBJFILES))
INCS = $(addprefix $(SRCDIR),$(INCLFILES))

CFLAGS = -fopenmp -std=c99 -fstrict-aliasing $(OPTIMIZE) $(GSL_INCL) $(MPI_INCL) $(FFT_INCL)

LIBS = -lm -lgsl -lgslcblas $(MPI_LIBS) $(GSL_LIBS) $(FFTW_LIBS)

# rules

%.o: %.c
	@echo [CC] $@
	@$(CC) $(CFLAGS)  -o $@ -c $<

$(EXEC)	: $(OBJS)
	@echo " "
	@echo 'SYSTYPE =' $(SYSTYPE)
	@echo 'CFLAGS =' $(CFLAGS)
	@echo 'LDFLAGS =' $(LIBS)
	@echo 'EXEC =' $(EXEC)
	@echo " "
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	@cd src && ctags -R *.[ch]

$(OBJS)	: $(INCS)

$(SRCDIR)config.h : Config 
	@echo [CC] $(SRCDIR)config.h
	@sed '/^#/d; /^$$/d; s/^/#define /g' Config > $(SRCDIR)config.h

$(SRCDIR)print_settings.c : Config	# does not work with sh shell
	@echo [CC] $(SRCDIR)print_settings.c
	@echo '/* Autogenerated File  */' >  $(SRCDIR)print_settings.c
	@echo '#include "globals.h"' >>  $(SRCDIR)print_settings.c
	@echo '#include "proto.h"' >>  $(SRCDIR)print_settings.c
	@echo 'void Print_compile_time_settings(){' >> $(SRCDIR)print_settings.c
	@echo '	rprintf("Compiled with : \n"' >> $(SRCDIR)print_settings.c
	@sed '/^#/d; /^$$/d; s/^/"      /g; s/$$/ \\n"/g;' Config >> \
		$(SRCDIR)print_settings.c
	@echo '); return ;}' >> $(SRCDIR)print_settings.c

clean : 
	rm -f $(OBJS) $(EXEC) src/config.h src/print_settings.c
