b:=../build
modflags := -I$(b) -J$(b)
flags    := -Wall
fc := gfortran
compiler := $(fc) $(modflags) $(flags)
#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)
lb := build

print-o := $(lb)/printer_mod.o
lib-a   := $(b)/libscme.a
gen-f90 := generic_program.f90
detr-o  := $(b)/detrace_mod.o


vpath %.f90 ../src

all:
	make -j4 run.x
	make -j4 tests.x

rrun:
	make -j4 run.x && ./run.x

ttests:
	make -j4 tests.x && ./tests.x
	

prog.x: poldip.f90 $(lib-a)
	$(compiler) $^ -o $@

stupid.x:$(detr-o) $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -D'STUPID'

detraceo.x: $(detr-o) $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) -DSUB=io_octa


detraceh.x: $(detr-o) $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) -DSUB=io_hexa


#$(detr-o) $(b)/detrace_mod.o
detraceg.x: $(detr-o) $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) -DSUB=io_general

apple.x: $(b)/detrace_apple.o $(gen-f90) 
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

qpole.x: $(b)/qpole_dev.o $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

opole.x: $(b)/opole_dev.o $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))


good := $(lb)/compressed_arrays.o generic_program.f90 $(lb)/printer_mod.o

test_good := $(good) $(lb)/compressed_foolin.o $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o $(lb)/data_types.o 


utils.x: $(lb)/compressed_utils.o  $(good)
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

convert.x: $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(good)
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

fool.x: $(lb)/compressed_foolin.o $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(good) 
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))


tests.x: $(lb)/compressed_tests.o $(test_good)
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))


run.x: $(lb)/compressed_run.o $(lb)/compressed_tests.o $(test_good)
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))



$(lb)/%.o:%.f90
	gfortran -I$(lb) -J$(lb) -c $< -o $@ -Wall

$(b)/%.o:%.f90
	$(compiler) -c $^ -o $@


$(lib-a):$(wildcard ../src/*.f90)
	(cd .. && make)

# Hierarchical dependencies

$(lb)/compressed_utils.o: $(print-o) $(lb)/compressed_arrays.o

$(lb)/compressed_tensors.o: $(lb)/compressed_utils.o 

$(lb)/compressed_foolin.o: $(lb)/compressed_tensors.o 

$(lb)/compressed_tests.o: $(lb)/compressed_foolin.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o 

$(lb)/compressed_run.o: $(lb)/compressed_tests.o 
# compressed run actually depend on all of the above objects, but if A list B, and B lists C, no need for A to list C

$(lb)/calc_derivs.o:$(lb)/sf_disp_tangtoe.o
$(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o: $(lb)/data_types.o


clean:
	rm *.x build/*

comp-clean:
	rm *.x build/comp*

