########################################################################
# New make file by JÖ. 
# This is a quite simple makefile that first makes sure that the modules are made by serially checking/creating them in the listed order. 
# Then the objects are comiled in parallel independent of order. 

# gfortran doesn't update good module files even if the corresponding object file is updated, therfore "touch $@" below. 
# otherwise the old module file causes this command to execute every time. 
# the alternative would be to make modules separately at from-scratch compiles, and thus separate the "makemodules" 
# and "makeobjects" into two different commands. 
########################################################################

# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:
.PHONY: clean modules-only objects-only all


######################################### directories:  
BD = build
SRCD = src

######################################### source paths:  
vpath %.f90 $(SRCD)
vpath %.cpp $(SRCD)


######################################### Compilers:  
FC = gfortran
CC = g++
######################################### Compiler flags:  
FFLAGS = -I$(BD) -J$(BD)   -pg -fPIC           -Ofast -march=native  -msoft-float -mavx  # 
CFLAGS = -I$(BD) -J$(BD)   -lstdc++            -Ofast -march=native #-O2 
#### module directory^ ### important^ ######### optimization^ (slower compile, faster code)

######################################### filenames in order:  
files_by_dependency = ${addprefix $(BD)/,  \
	data_types \
	parameters \
	max_parameters \
	multipole_parameters \
	polariz_parameters \
	calcEnergy_mod \
	tang_toennies \
	molecProperties \
	calc_derivs \
	calc_higher_order \
	calc_lower_order \
	inducePoles \
	forceCM_mod \
	torqueCM_mod \
	mdutil \
	molforce \
	atomicForces_mod \
	dispersion_mod \
	rho \
	coreInt_mod \
	scme   }

######################################### path/filename.ext :
fort_obj = $(files_by_dependency:=.o)
fort_mod = $(files_by_dependency:=.mod)
cpp_obj  = $(BD)/ps.o

######################################### easy commands:
all:
	make -j1 modules-only
	make -j4 objects-only

clean:
	rm $(addprefix $(BD)/, *.o *.a *.mod)


######################################### making modules:
modules-only:$(fort_mod)


$(fort_mod):$(BD)/%.mod: %.f90
	$(FC) $(FFLAGS) -fsyntax-only $<; touch $@
	

######################################### making objects:
objects-only: $(BD)/libscme.a

$(BD)/libscme.a: $(cpp_obj) $(fort_obj) $(fort_mod)
	ar rcs $@ $^ 

$(fort_obj):$(BD)/%.o: %.f90 
	$(FC) $(FFLAGS) -c $< -o $@

$(cpp_obj):$(BD)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

