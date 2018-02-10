### DIRECTORIES
b:=../build
ob := obuild
lb := lbuild
dodir=@mkdir -p $(@D)


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



### COMPILATION 
fc := gfortran
opti = 
#-Ofast -flto -march=native 
# -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -finline-functions 
err := -Wall -fmax-errors=3 -fcheck=all
#-Wall -fmax-errors=3 -fcheck=all
# -fcheck=all
# $(err)
lib_comp := $(fc) -I$(b) -J$(ob) $(err) $(opti) -ffree-line-length-0
loc_comp := $(fc) -I$(lb) -J$(lb) $(err) $(opti) -ffree-line-length-0
# Must place main object first for generic program to call the correct module. 
generic = $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) generic_program.f90
comp1 = $(lib_comp) $(generic)
comp2 = $(loc_comp) $(generic)

vpath %.f90 ../src

# -I = where to search! (for #/include and .mod)
# -J = where to put!    (.mod)


$(lb)/%.o:%.f90
	$(dodir)
	$(loc_comp) -c $< -o $@

$(ob)/%.o:%.f90
	$(dodir)
	$(fc) -J$(ob) -I$(ob) -I$(b) -c $^ -o $@  -Wall

$(lb)/detrace_apple.o:detrace_apple.f90
	$(dodir)
	$(fc) -I$(lb) -J$(lb) -c $< -o $@


$(lib-a):$(wildcard ../src/*.f90)
	(cd .. && make)





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
	rm *.x lbuild/comp*

