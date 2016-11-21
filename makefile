# New make file by JÖ. 
# The previous one dealt with "compiling modules before their depdendencies" through stating dependencies between object files. This resulted in a complicated make-file and longer compile-times for small changes since unecessayily large parts of the program were recomilesd. Now this issue is solved through ordering the file names instead. This might seem less robust. It might be. 


# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:

BD = build
SRCD = src

vpath %.f90 $(SRCD)
vpath %.cpp $(SRCD)
#vpath %.mod $(BD)


FC = gfortran
CC = g++
FFLAGS = -pg  -I$(BD) -J$(BD) -fPIC -Ofast -march=native  -msoft-float -mavx  #
CFLAGS = -O2 -I$(BD) -J$(BD) -lstdc++ -Ofast -march=native
#-O3

######################################### Files in order:  
fort_obj = $(addprefix $(BD)/,  \
	data_types.o parameters.o max_parameters.o \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o tang_toennies.o molecProperties.o calc_higher_order.o calc_lower_order.o \
	inducePoles.o calc_derivs.o forceCM_mod.o torqueCM_mod.o \
	mdutil.o molforce.o atomicForces_mod.o \
	dispersion_mod.o rho.o coreInt_mod.o \
	scme.o )

cpp_obj = $(BD)/ps.o

######################################### Compile:

all: $(BD)/libscme.a

$(BD)/libscme.a: $(fort_obj) $(cpp_obj)
	ar rcs $@ $^

$(BD)/%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $<

$(BD)/%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<



######################################### Clean:
.PHONY: clean
clean:
	rm $(addprefix $(BD)/, *.o *.a *.mod)



######################################### Trash:
#$(BD)/%.o:$(BD)/%.mod
#$(BD)/%.mod: %.f90
#	$(FC) $(FFLAGS) -fsyntax-only $<
