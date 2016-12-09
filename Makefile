# **********************************************************************************************************************************
# settings

# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:

OBJDIR = obj
MODDIR = mod
SRCDIR = src/

vpath %.f90 $(SRCDIR)
vpath %.cpp $(SRCDIR)

FC = gfortran
CC = g++
opti = -O2
# -Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native 

FFLAGS = $(opti) -pg -I$(MODDIR) -J$(MODDIR)
CFLAGS = -O2 -I$(MODDIR) -J$(MODDIR) -lstdc++
#-fopenmp
#-floop-unroll-and-jam -ftree-loop-if-convert
vect = -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution

OBJ = $(addprefix $(OBJDIR)/, \
	scme.o calc_derivs.o calc_higher_order.o \
	data_types.o parameters.o max_parameters.o \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o calc_lower_order.o \
	inducePoles.o forceCM_mod.o torqueCM_mod.o \
	atomicForces_mod.o molforce.o tang_toennies.o mdutil.o \
	molecProperties.o dispersion_mod.o coreInt_mod.o rho.o ps.o)
#OBJC = $(addprefix $(OBJDIR)/, ps.o)
HEADERS = $(addprefix $(OBJDIR)/, constants.h ps.h)

#all: $(OBJDIR)/scme.o
all:
	make -j4 it

it:$(OBJDIR)/libscme.a

# **********************************************************************************************************************************

# linking
#

# library
$(OBJDIR)/libscme.a: $(OBJ)
	ar rcs $@ $^

# compiling

#$(OBJDIR)/.f90.o:
#	$(FC) $(FFLAGS) -c -o $@ $<

#$(OBJDIR)/.cpp.o: $(HEADERS)
#	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

#clean:

######################################### Dependencies:
# Single depdenden --- multiple prerequisites:

$(OBJDIR)/scme.o:		\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/parameters.o		\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/forceCM_mod.o		\
$(OBJDIR)/torqueCM_mod.o	\
$(OBJDIR)/atomicForces_mod.o	\
$(OBJDIR)/calcEnergy_mod.o	\
$(OBJDIR)/coreInt_mod.o		\
$(OBJDIR)/dispersion_mod.o	\
$(OBJDIR)/multipole_parameters.o\

$(OBJDIR)/atomicForces_mod.o:	\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/molforce.o		\


$(OBJDIR)/molforce.o:	\
$(OBJDIR)/data_types.o	\
$(OBJDIR)/mdutil.o	\

$(OBJDIR)/coreInt_mod.o:	\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/parameters.o 		\
$(OBJDIR)/rho.o			\

$(OBJDIR)/rho.o:		\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/parameters.o		\

$(OBJDIR)/molecProperties.o:	\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o 	\
$(OBJDIR)/tang_toennies.o	\

$(OBJDIR)/tang_toennies.o:	\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/parameters.o		\



# multiple dependents --- few prerequisite:
$(OBJDIR)/mdutil.o		\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/multipole_parameters.o\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/parameters.o		\
:$(OBJDIR)/data_types.o		\


$(OBJDIR)/dispersion_mod.o	\
$(OBJDIR)/torqueCM_mod.o	\
$(OBJDIR)/forceCM_mod.o		\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/calcEnergy_mod.o	\
:$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\


$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
:$(OBJDIR)/data_types.o		\
$(OBJDIR)/max_parameters.o	\
$(OBJDIR)/molecProperties.o	\


# **********************************************************************************************************************************
# module dependencies

#$(OBJDIR)/parameters.o: $(OBJDIR)/data_types.o
#
#$(OBJDIR)/tang_toennies.o: $(OBJDIR)/parameters.o
#
#$(OBJDIR)/rho.o:$(OBJDIR)/parameters.o
#
#$(OBJDIR)/scme.o: $(addprefix $(OBJDIR)/, data_types.o parameters.o max_parameters.o multipole_parameters.o polariz_parameters.o \
#	molecProperties.o calc_higher_order.o calc_lower_order.o calc_derivs.o inducePoles.o forceCM_mod.o \
#	torqueCM_mod.o atomicForces_mod.o calcEnergy_mod.o coreInt_mod.o dispersion_mod.o ps.o)
#
#$(OBJDIR)/molecProperties.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/tang_toennies.o
#
#$(OBJDIR)/calc_higher_order.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/molecProperties.o
#
#$(OBJDIR)/calc_lower_order.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/molecProperties.o
#
#$(OBJDIR)/calc_derivs.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/molecProperties.o
#
#$(OBJDIR)/coreInt_mod.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/parameters.o $(OBJDIR)/rho.o
#
#$(OBJDIR)/atomicForces_mod.o: $(OBJDIR)/data_types.o $(OBJDIR)/max_parameters.o $(OBJDIR)/molforce.o
#
#$(OBJDIR)/molforce.o: $(OBJDIR)/data_types.o $(OBJDIR)/mdutil.o

# **********************************************************************************************************************************
# cleanup

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/*.a $(MODDIR)/*.mod


# **********************************************************************************************************************************
