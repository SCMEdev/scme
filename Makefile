# ******************************************************************************************************
# settings

# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:

OBJDIR = obj
MODDIR = mod
SRCDIR = src/
NEW = new

#have added the "new/" directory and PS-files
vpath %.f90 $(SRCDIR)
vpath %.f90 $(NEW)
vpath %.cpp $(SRCDIR)

FC = gfortran
CC = g++
FFLAGS = -O3 -I$(MODDIR) -J$(MODDIR)
CFLAGS = -O2 -I$(MODDIR) -J$(MODDIR) -lstdc++
## Flags for mor optimizations and 49 passes at home:
#-Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native -fopenmp -finline-functions
# -Ofast
FFLAGS = $(opti) -pg -I$(MODDIR) -J$(MODDIR)
CFLAGS = $(opti) -I$(MODDIR) -J$(MODDIR) -lstdc++
#-fopenmp
#-floop-unroll-and-jam -ftree-loop-if-convert
#vect = -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution

OBJ = $(addprefix $(OBJDIR)/, \
	scme_ps.o calc_derivs.o calc_higher_order.o \
	data_types.o parameters.o \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o calc_lower_order.o \
	inducePoles.o \
	atomicForces_mod.o molforce.o mdutil.o \
	molecProperties.o \
	ps_pes.o ps_dms.o constants.o printer_mod.o sf_disp_tangtoe.o force_torqueCM.o)
	
#	 max_parameters.o forceCM_mod.o torqueCM_mod.o tang_toennies.o rho.o dispersion_mod.o coreInt_mod.o ps.o
#OBJC = $(addprefix $(OBJDIR)/, ps.o)
#HEADERS = $(addprefix $(OBJDIR)/, constants.h ps.h)

#all: $(OBJDIR)/scme.o
all:
	make -j4 it

it:$(OBJDIR)/libscme.a

# *****************************************************************************************************************

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

#$(OBJDIR)/%.o: %.cpp
#	$(CC) $(CFLAGS) -c -o $@ $<

#/////////////////////////////////////////////////// Dependencies //////
# special dependencies:
$(OBJDIR)/sf_disp_tangtoe.o:$(OBJDIR)/parameters.o
$(OBJDIR)/atomicForces_mod.o:$(OBJDIR)/molforce.o
$(OBJDIR)/molforce.o:$(OBJDIR)/mdutil.o
$(OBJDIR)/ps_dms.o	$(OBJDIR)/ps_pes.o:$(OBJDIR)/constants.o

$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o:	\
$(OBJDIR)/sf_disp_tangtoe.o \


# scme dep. on most
$(OBJDIR)/scme_ps.o:		\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/parameters.o		\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/atomicForces_mod.o	\
$(OBJDIR)/calcEnergy_mod.o	\
$(OBJDIR)/multipole_parameters.o\
$(OBJDIR)/ps_pes.o \
$(OBJDIR)/ps_dms.o \
$(OBJDIR)/printer_mod.o \
$(OBJDIR)/sf_disp_tangtoe.o \
$(OBJDIR)/force_torqueCM.o \

# most dep. on data_types
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
$(OBJDIR)/printer_mod.o \
$(OBJDIR)/force_torqueCM.o \
$(OBJDIR)/sf_disp_tangtoe.o \
$(OBJDIR)/ps_dms.o	\
$(OBJDIR)/ps_pes.o \
$(OBJDIR)/molforce.o	\
$(OBJDIR)/atomicForces_mod.o	\
$(OBJDIR)/force_torqueCM.o		\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/calcEnergy_mod.o	\
$(OBJDIR)/mdutil.o		\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/multipole_parameters.o\
$(OBJDIR)/parameters.o:		\
$(OBJDIR)/data_types.o		\



.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(OBJDIR)/*.a $(MODDIR)/*.mod


# *********************************************************************************************************************
