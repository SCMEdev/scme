# New make file by JÖ. 
# The previous one dealt with "compiling modules before their depdendencies" through stating dependencies between object files. 
# This resulted in a complicated make-file and longer compile-times for small changes since unecessayily large parts of the program were recomilesd. 
# Now this issue is solved through ordering the file names instead. This might seem less robust. It might be. 


# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:
#.SUFFIXES: .f90 .o .cpp

BD = build
SRCD = src

vpath %.f90 $(SRCD)
vpath %.cpp $(SRCD)
#vpath %.mod $(BD)


FC = gfortran
CC = g++
FFLAGS = -I$(BD) -J$(BD)  #-pg -fPIC -Ofast -march=native  -msoft-float -mavx  #
CFLAGS = -I$(BD) -J$(BD) #-O2 -lstdc++ -Ofast -march=native
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

######################################### 

all:
	make -j1 makemodules
	make -j4 makeobjects

makemodules:$(BD)/scme.mod

$(BD)/%.mod: %.f90
	$(FC) $(FFLAGS) -fsyntax-only $<; touch $@
	


makeobjects: $(BD)/libscme.a

$(BD)/libscme.a: $(fort_obj) $(cpp_obj) 
	ar rcs $@ $^ 


$(BD)/%.o: %.f90 
	$(FC) $(FFLAGS) -c $< -o $@

$(BD)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

#$(BD)/%.mod: %.f90
#	$(FC) $(FFLAGS) -fsyntax-only $< 




######################################### Compile:

$(BD)/scme.mod\
:$(BD)/data_types.mod\
$(BD)/max_parameters.mod\
$(BD)/parameters.mod\
$(BD)/polariz_parameters.mod\
$(BD)/molecProperties.mod\
$(BD)/calc_lower_order.mod\
$(BD)/calc_higher_order.mod\
$(BD)/calc_derivs.mod\
$(BD)/inducePoles.mod\
$(BD)/forceCM_mod.mod\
$(BD)/torqueCM_mod.mod\
$(BD)/atomicForces_mod.mod\
$(BD)/calcEnergy_mod.mod\
$(BD)/coreInt_mod.mod\
$(BD)/dispersion_mod.mod\
$(BD)/multipole_parameters.mod

$(BD)/atomicForces_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/molforce.mod	


$(BD)/molforce.mod:$(BD)/data_types.mod $(BD)/mdutil.mod 


$(BD)/mdutil.mod:$(BD)/data_types.mod 			
$(BD)/polariz_parameters.mod:$(BD)/data_types.mod 			
$(BD)/multipole_parameters.mod:$(BD)/data_types.mod 			
$(BD)/max_parameters.mod:$(BD)/data_types.mod 			
$(BD)/parameters.mod:$(BD)/data_types.mod 			



$(BD)/dispersion_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 
$(BD)/torqueCM_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 
$(BD)/forceCM_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 
$(BD)/inducePoles.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 
$(BD)/tang_toennies.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 
$(BD)/calcEnergy_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod 

$(BD)/coreInt_mod.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/parameters.mod $(BD)/rho.mod
$(BD)/rho.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/parameters.mod

$(BD)/molecProperties.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/tang_toennies.mod 

$(BD)/calc_derivs.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/molecProperties.mod 
$(BD)/calc_lower_order.mod:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/molecProperties.mod 


######################################### Clean:
.PHONY: clean
clean:
	rm $(addprefix $(BD)/, *.o *.a *.mod)

remake: 
	make clean
	make

######################################### 

#$(BD)/mdutil.mod 				\
#$(BD)/polariz_parameters.mod 	\
#$(BD)/multipole_parameters.mod 	\
#$(BD)/max_parameters.mod 		\
#$(BD)/parameters.mod			\
#:$(BD)/data_types.mod 			
#
#$(BD)/atomicForces_mod.mod 		\
#:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/molforce.mod	
#
#
#$(BD)/molforce.mod \
#:$(BD)/data_types.mod $(BD)/mdutil.mod 
#
#
#$(BD)/dispersion_mod.mod	\
#$(BD)/torqueCM_mod.mod		\
#$(BD)/forceCM_mod.mod 		\
#$(BD)/inducePoles.mod 		\
#$(BD)/tang_toennies.mod 	\
#$(BD)/calcEnergy_mod.mod 	\
#:$(BD)/data_types.mod $(BD)/max_parameters.mod 
#
#$(BD)/coreInt_mod.mod	\
#$(BD)/rho.mod 			\
#:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/parameters.mod
#
#$(BD)/molecProperties.mod \
#:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/tang_toennies.mod 
#
#$(BD)/calc_derivs.mod 		\
#$(BD)/calc_lower_order.mod 	\
#$(BD)/calc_higher_order.mod \
#:$(BD)/data_types.mod $(BD)/max_parameters.mod $(BD)/molecProperties.mod 
#
#
#$(BD)/scme.mod \
#:$(BD)/data_types.mod 			\
#$(BD)/parameters.mod 			\
#$(BD)/max_parameters.mod 		\
#$(BD)/multipole_parameters.mod	\
#$(BD)/polariz_parameters.mod 	\
#$(BD)/calcEnergy_mod.mod 		\
#$(BD)/tang_toennies.mod 		\
#$(BD)/molecProperties.mod 		\
#$(BD)/calc_higher_order.mod 	\
#$(BD)/calc_lower_order.mod 		\
#$(BD)/inducePoles.mod 			\
#$(BD)/calc_derivs.mod 			\
#$(BD)/forceCM_mod.mod 			\
#$(BD)/torqueCM_mod.mod			\
#$(BD)/mdutil.mod 				\
#$(BD)/molforce.mod 				\
#$(BD)/atomicForces_mod.mod 		\
#$(BD)/dispersion_mod.mod 		\
#$(BD)/rho.mod 					\
#$(BD)/coreInt_mod.mod	


