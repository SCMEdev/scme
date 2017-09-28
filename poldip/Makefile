scme_build=../build
modflags = -I$(scme_build) -J$(scme_build)
flags = -Wall

lib = $(scme_build)/libscme.a
#$(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o qpole.o opole.o polariz_parameters.o)

prog.x: get_poldip.f90 $(lib)
	gfortran $(modflags) $(flags) $^ -o $@

#	gfortran ~/gits/scme/src/multipole_parameters.f90 ~/gits/scme/src/localAxes_mod.f90 ~/gits/scme/src/data_types.f90 ~/gits/scme/src/printer_mod.f90 qpole.f90

#$(lib):
#	(cd ../ && make)

clean:
	rm prog.x

$(lib):$(wildcard ../src/*.f90)
	(cd .. && make)
