# disable the built-in (implicit) rules to avoid trying to compile X.o from X.mod (Modula-2 program)
.SUFFIXES:

OBJDIR = obj
MODDIR = mod
SRCDIR = src/
#NEW = new

dirs = $(OBJDIR) $(MODDIR)

#have added the "new/" directory and PS-files
vpath %.f90 $(SRCDIR)
#vpath %.f90 $(NEW)
#vpath %.cpp $(SRCDIR)

FC = gfortran
CC = g++
opti = -O0

## optimization:
#opti = -Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native -fopenmp -finline-functions
## warn all:
#-Wall

FFLAGS = $(opti) -pg -I$(MODDIR) -J$(MODDIR) 
CFLAGS = $(opti) -I$(MODDIR) -J$(MODDIR) -lstdc++


OBJ = $(addprefix $(OBJDIR)/, \
	scme_ps.o calc_derivs.o calc_higher_order.o \
	data_types.o \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o calc_lower_order.o \
	inducePoles.o \
	molecProperties.o \
	ps_pes.o ps_dms.o printer_mod.o sf_disp_tangtoe.o force_torqueCM.o \
	localAxes_mod.o )
	
#	 max_parameters.o forceCM_mod.o torqueCM_mod.o tang_toennies.o rho.o dispersion_mod.o coreInt_mod.o ps.o parameters.o constants.o molforce.o mdutil.o 	atomicForces_mod.o \

#OBJC = $(addprefix $(OBJDIR)/, ps.o)
#HEADERS = $(addprefix $(OBJDIR)/, constants.h ps.h)

#/// Build

all:
	make -j4 it

it:$(OBJDIR)/libscme.a

# library

$(OBJDIR)/libscme.a: $(OBJ) $(dirs)
	ar rcs $@ $(OBJ)

$(dirs):
	mkdir $@


$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

#//// Clean
.PHONY: clean
clean:
	rm -f $(OBJDIR)/* $(MODDIR)/*


#/// Dependencies

# special dependencies:
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o:	\
$(OBJDIR)/sf_disp_tangtoe.o \


# scme dep. on most
$(OBJDIR)/scme_ps.o:		\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/data_types.o		\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/calcEnergy_mod.o	\
$(OBJDIR)/multipole_parameters.o\
$(OBJDIR)/ps_pes.o \
$(OBJDIR)/ps_dms.o \
$(OBJDIR)/printer_mod.o \
$(OBJDIR)/sf_disp_tangtoe.o \
$(OBJDIR)/force_torqueCM.o \
$(OBJDIR)/localAxes_mod.o \



$(OBJDIR)/localAxes_mod.o \
$(OBJDIR)/molecProperties.o	\
$(OBJDIR)/calc_derivs.o		\
$(OBJDIR)/calc_lower_order.o	\
$(OBJDIR)/calc_higher_order.o	\
$(OBJDIR)/printer_mod.o \
$(OBJDIR)/force_torqueCM.o \
$(OBJDIR)/sf_disp_tangtoe.o \
$(OBJDIR)/ps_dms.o	\
$(OBJDIR)/ps_pes.o \
$(OBJDIR)/force_torqueCM.o		\
$(OBJDIR)/inducePoles.o		\
$(OBJDIR)/calcEnergy_mod.o	\
$(OBJDIR)/polariz_parameters.o	\
$(OBJDIR)/multipole_parameters.o: \
$(OBJDIR)/data_types.o		\




