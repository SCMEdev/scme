b:=../build
modflags := -I$(b) -J$(b)
flags    := -Wall
compiler := gfortran $(modflags) $(flags)
#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)


lib-a   := $(b)/libscme.a
gen-f90 := generic_program.f90
detr-o  := $(b)/detrace_mod.o


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

convert.x: $(b)/compressed_tensors.o $(gen-f90) $(lib-a)
	$(compiler) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

#-cpp -D'STUPID'

$(b)/%.o:%.f90
	$(compiler) -c $^ -o $@
# -cpp

$(lib-a):$(wildcard ../src/*.f90)
	(cd .. && make)


clean:
	rm -f *.x *.o *.mod

