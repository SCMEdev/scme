#.SUFIXES:
#.PREFIXES

### VARIABLES
b:=build
lib:=$b/libscme.a

vpath %.f90 src
vpath %.f90 devdir

FC:=gfortran
#FC:=flang
opti:= -O0
#opti = -Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native -finline-functions -fopenmp 

FFLAGS:= $(opti) -pg -I$b -J$b -cpp -D'DEBUG_PRINTING' -ffree-line-length-0 -Wall
#-Wall
#CFLAGS:= $(opti) -I$b -J$b -lstdc++

### OBJECT FILES grouped by dependencies
# Depend only on data_types:
bulk_obj:=$(addprefix $b/, \
	multipole_parameters.o polariz_parameters.o \
	calcEnergy_mod.o inducePoles.o ps_pes.o ps_dms.o printer_mod.o \
	sf_disp_tangtoe.o force_torqueCM.o localAxes_mod.o qpole.o opole.o )

# Depend also on swiching-func (sf):
sf_obj:= $(addprefix $b/, \
	calc_derivs.o calc_higher_order.o calc_lower_order.o molecProperties.o )
	
# Depends on all/none (ENDpoints of dep.-tree)
end_obj:= $(addprefix $b/, \
	scme_ps.o data_types.o )

# extra from development dir
dev_obj := $(addprefix $b/, \
	compressed_utils.o detrace_apple.o compressed_arrays.o )
	# compressed_tensors.o

### OBJECT LISTS:
most_obj:=$(sf_obj) $(bulk_obj)
all_obj:= $(dev_obj) $(most_obj) $(end_obj)



### RULES:
all:
	make -r -j4 $(lib)

# (.mod-files still needed for compilation)
$(lib): $(all_obj) 
	ar rcs $@ $^


$b/%.o: %.f90 
	$(FC) $(FFLAGS) -o $@ -c $<

$b:
	mkdir $@

.PHONY: clean
clean:
	rm -f $b/*


### DEPENDENCIES

#temporary >>>
$b/compressed_utils.o:$b/compressed_arrays.o $b/printer_mod.o
$b/calc_derivs.o : $b/detrace_apple.o
$b/scme_ps.o:  $b/compressed_utils.o
#temporary <<<

$(sf_obj): $b/sf_disp_tangtoe.o 
$(most_obj): $b/data_types.o
$b/scme_ps.o:  $(most_obj) 
