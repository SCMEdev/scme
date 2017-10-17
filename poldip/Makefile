scme_build=../build
modflags = -I$(scme_build) -J$(scme_build)
flags = -Wall
compiler = gfortran $(modflags) $(flags)
lib = $(scme_build)/libscme.a
#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)

prog.x: poldip.f90 $(lib)
	$(compiler) $^ -o $@

test_det.x: detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'RUN_TEST'

other_det.x: detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'RUN_OTHER'

detraceh.x: detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'HEXAD'

detraceo.x: detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'OCTA'

detraceg.x: detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'GENERAL'

stupid.x:detrace_mod.f90 detrace_program.f90 $(lib)
	$(compiler) $^ -o $@ -cpp -D'STUPID'


#	gfortran ~/gits/scme/src/multipole_parameters.f90 ~/gits/scme/src/localAxes_mod.f90 ~/gits/scme/src/data_types.f90 ~/gits/scme/src/printer_mod.f90 qpole.f90

#$(lib):
#	(cd ../ && make)

clean-prog:
	rm prog.x

clean-det:
	rm detrace.x

clean-oth:
	rm other_det.x

clean-progs:
	rm *.x

$(lib):$(wildcard ../src/*.f90)
	(cd .. && make)
