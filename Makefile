### VARIABLES
b:=build
s:=src

FC:=gfortran
#FC:=flang
opti:= -O0
#opti = -Ofast -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -march=native -fopenmp -finline-functions
#-Wall

FFLAGS:= $(opti) -pg -I$b -J$b -cpp -D'DEBUG_PRINTING' 
#CFLAGS:= $(opti) -I$b -J$b -lstdc++

### FILE NAMES (according to dependencies)
# Depend only on data_types:
bulk_names:=\
	multipole_parameters polariz_parameters \
	calcEnergy_mod inducePoles ps_pes ps_dms printer_mod \
	sf_disp_tangtoe force_torqueCM localAxes_mod qpole opole\
	detrace_apple

# Depend also on swiching-func (sf):
sf_names:=\
	calc_derivs calc_higher_order\
	calc_lower_order molecProperties

print_names :=\
	compressed_tensors

# Depends on all/none (ENDpoints of dep.-tree)
end_names:=\
	scme_ps data_types

### OBJECT LISTS:
sf_obj   := $(addprefix $b/, $(sf_names:=.o) )
print_obj:= $(addprefix $b/, $(print_names:=.o) )
most_obj := $(sf_obj) $(print_obj)  $(addprefix $b/, $(bulk_names:=.o))
all_obj  := $(most_obj) $(addprefix $b/, $(end_names:=.o) ) 

lib:=$b/libscme.a


### RULES:
all:
	make -j4 $(lib)

# (.mod-files still needed for compilation)
$(lib): $(all_obj) 
	ar rcs $@ $^

#temporary !!!!!!!!!/////////////////////////////////////////////////////////////////////
$b/compressed_tensors.o:devdir/compressed_tensors.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$b/detrace_apple.o:devdir/detrace_apple.f90
	$(FC) $(FFLAGS) -o $@ -c $<

$b/%.o: $s/%.f90 
	$(FC) $(FFLAGS) -o $@ -c $<

$b:
	mkdir $@

.PHONY: clean
clean:
	rm -f $b/*


### DEPENDENCIES
# swich-func dep:
$(sf_obj): $b/sf_disp_tangtoe.o 

$(print_obj): $b/printer_mod.o


$/scme_ps.o $b/calc_derivs.o: $b/detrace_apple.o $b/compressed_tensors.o
# scme depends on all:
$b/scme_ps.o: $(most_obj)

# all depend on data_types:
$(most_obj) $b/scme_ps.o: $b/data_types.o

# printer dep:

