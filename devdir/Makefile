scme_build:=../build
modflags := -I$(scme_build) -J$(scme_build)
flags := -Wall
compiler := gfortran $(modflags) $(flags)
lib := $(scme_build)/libscme.a
#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)

genericp := generic_program.f90

obj := detrace_mod.o


prog.x: poldip.f90 $(lib)
	$(compiler) $^ -o $@

test_det.x: $(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'RUN_TEST'

other_det.x: $(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'RUN_OTHER'

detraceh.x: $(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'HEXAD'

detraceo.x: $(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'OCTA'

detraceg.x: $(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'GENERAL'

stupid.x:$(obj) $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -D'STUPID'

apple.x: detrace_apple.o $(genericp) 
	$(compiler) $^ -o $@ -cpp -DAPPLE


qpole.x: qpole_dev.o $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -DYO=$(<:.o=)

opole.x: opole_dev.o $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -DYO=$(<:.o=)

symmetrize.x: symmetrize.o $(genericp) $(lib)
	$(compiler) $^ -o $@ -cpp -DYO=$(<:.o=)

#-cpp -D'STUPID'

%.o:%.f90
	$(compiler) -c $^ -o $@ -cpp

$(lib):$(wildcard ../src/*.f90)
	(cd .. && make)


clean:
	rm -f *.x *.o *.mod

