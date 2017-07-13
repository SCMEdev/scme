scme_build=../build
modflags = -I$(scme_build) -J$(scme_build)

objs = $(addprefix $(scme_build)/, multipole_parameters.o localAxes_mod.o data_types.o printer_mod.o)

prog.x: $(objs) qpole.f90
	gfortran $(modflags) $^ -o $@

#	gfortran ~/gits/scme/src/multipole_parameters.f90 ~/gits/scme/src/localAxes_mod.f90 ~/gits/scme/src/data_types.f90 ~/gits/scme/src/printer_mod.f90 qpole.f90

$(objs):
	(cd ../ && make)
