### DIRECTORIES
b:=../build
ob := obuild
lb := lbuild


### COMPILATION 
fc := gfortran
lib_comp := $(fc) -I$(b) -J$(ob)
loc_comp := $(fc) -I$(lb) -J$(lb)
# Must place main object first for generic program to call the correct module. 
generic = $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) generic_program.f90
comp1 = $(lib_comp) $(generic)
comp2 = $(loc_comp) $(generic)

$(lb)/%.o:%.f90
	$(loc_comp) -c $< -o $@

$(ob)/%.o:%.f90
	$(fc) -J$(ob) -I$(ob) -I$(b) -c $^ -o $@  -Wall


$(lib-a):$(wildcard ../src/*.f90)
	(cd .. && make)

vpath %.f90 ../src

# -I = where to search! (for #/include and .mod)
# -J = where to put!    (.mod)






# Fast commands
all:
	make -j4 run.x
	make -j4 tests.x

rrun:
	make -j4 run.x && ./run.x

ttests:
	make -j4 tests.x && ./tests.x

old: 
	make -j4 prog.x stupid.x detraceo.x detraceh.x detraceg.x qpole.x opole.x

new: 
	make -j4 apple.x utils.x convert.x fool.x tests.x run.x 





### OLDER

print-o := $(ob)/printer_mod.o
detr-o  := $(ob)/detrace_mod.o
lib-a   := $(b)/libscme.a

prog.x: $(ob)/poldip.o $(lib-a)
	$(lib_comp) $^ -o $@

stupid.x: generic_program.f90 $(lib-a)
	$(lib_comp) $^ -o $@ -cpp -D'STUPID'


detraceo.x: $(detr-o) $(lib-a)
	$(comp1) -DSUB=io_octa


detraceh.x: $(detr-o) $(lib-a)
	$(comp1) -DSUB=io_hexa


detraceg.x: $(detr-o) $(lib-a)
	$(comp1) -DSUB=io_general

qpole.x: $(ob)/qpole_dev.o $(lib-a)
	$(comp1)

opole.x: $(ob)/opole_dev.o $(lib-a)
	$(comp1)



### COMPRESSED

apple.x: $(lb)/detrace_apple.o 
	$(comp2)

objects := $(lb)/compressed_arrays.o $(lb)/printer_mod.o
# using 'objects':
utils.x: $(lb)/compressed_utils.o  $(objects)
	$(comp2)
convert.x: $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(objects)
	$(comp2)
fool.x: $(lb)/compressed_foolin.o $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(objects) 
	$(comp2)

test_objects := $(objects) $(lb)/compressed_foolin.o $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o $(lb)/data_types.o 
# using 'test_objects':
tests.x: $(lb)/compressed_tests.o $(test_objects)
	$(comp2)

run.x: $(lb)/compressed_run.o $(lb)/compressed_tests.o $(test_objects)
	$(comp2)


### RULES

### DEPENDENCIES

$(lb)/compressed_utils.o: $(lb)/printer_mod.o $(lb)/compressed_arrays.o

$(lb)/compressed_tensors.o: $(lb)/compressed_utils.o 

$(lb)/compressed_foolin.o: $(lb)/compressed_tensors.o 

$(lb)/compressed_tests.o: $(lb)/compressed_tensors.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o 

$(lb)/compressed_run.o: $(lb)/compressed_tests.o 

$(lb)/calc_derivs.o:$(lb)/sf_disp_tangtoe.o
$(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o: $(lb)/data_types.o


### CLEAN
clean:
	rm -f *.x *.mod *.o lbuild/* obuild/*

comp-clean:
	rm *.x build/comp*



#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)
