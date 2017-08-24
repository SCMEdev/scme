# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:

b = build
s = src
#NEW = new

#dirs = $(OBJDIR) $(MODDIR)
dirs = $b

#have added the "new/" directory and PS-files
vpath %.f90 $(s)
#vpath %.f90 $(NEW)
#vpath %.cpp $(SRCDIR)

FC = gfortran
CC = g++
opti = -O0

## optimization:
#opti = -Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native -fopenmp -finline-functions
## warn all:
#-Wall

#Debug prerpcessor flag only in this makefile, makes the program print stuff with the print routine:
FFLAGS = $(opti) -pg -I$b -J$b -cpp -D'DEBUG'
CFLAGS = $(opti) -I$b -J$b -lstdc++


OBJ = $(addprefix $b/, \
	scme_ps.o calc_derivs.o calc_higher_order.o \
	data_types.o \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o calc_lower_order.o \
	inducePoles.o \
	molecProperties.o \
	ps_pes.o ps_dms.o printer_mod.o sf_disp_tangtoe.o force_torqueCM.o \
	localAxes_mod.o qpole.o)
	
#	 max_parameters.o forceCM_mod.o torqueCM_mod.o tang_toennies.o rho.o dispersion_mod.o coreInt_mod.o ps.o parameters.o constants.o molforce.o mdutil.o 	atomicForces_mod.o \

#OBJC = $(addprefix $(OBJDIR)/, ps.o)
#HEADERS = $(addprefix $(OBJDIR)/, constants.h ps.h)

#/// Build

all:
	make -j4 it

it:$b/libscme.a

# library

$b/libscme.a: $(OBJ) $(dirs)
	ar rcs $@ $(OBJ) $b/*.mod

$b:
	mkdir $@


$b/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

#//// Clean
.PHONY: clean
clean:
	rm -f $b/*


#/// Dependencies

$b/qpole.o\
$b/localAxes_mod.o:\
$b/printer_mod.o 


# special dependencies:
$b/molecProperties.o	\
$b/calc_derivs.o		\
$b/calc_lower_order.o	\
$b/calc_higher_order.o:	\
$b/sf_disp_tangtoe.o \


# scme dep. on most
$b/scme_ps.o:		\
$b/calc_derivs.o		\
$b/data_types.o		\
$b/polariz_parameters.o	\
$b/molecProperties.o	\
$b/calc_lower_order.o	\
$b/calc_higher_order.o	\
$b/inducePoles.o		\
$b/calcEnergy_mod.o	\
$b/multipole_parameters.o\
$b/ps_pes.o \
$b/ps_dms.o \
$b/printer_mod.o \
$b/sf_disp_tangtoe.o \
$b/force_torqueCM.o \
$b/localAxes_mod.o \
$s/debug.h\
$b/qpole.o\

$b/qpole.o\
$b/localAxes_mod.o \
$b/molecProperties.o	\
$b/calc_derivs.o		\
$b/calc_lower_order.o	\
$b/calc_higher_order.o	\
$b/printer_mod.o \
$b/force_torqueCM.o \
$b/sf_disp_tangtoe.o \
$b/ps_dms.o	\
$b/ps_pes.o \
$b/force_torqueCM.o		\
$b/inducePoles.o		\
$b/calcEnergy_mod.o	\
$b/polariz_parameters.o	\
$b/multipole_parameters.o: \
$b/data_types.o		\

