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



convert.x: $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(lb)/detrace_apple.o $(print-o) $(gen-f90) $(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o $(lb)/data_types.o $(lb)/compressed_arrays.o
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))

utils.x: $(lb)/compressed_utils.o $(print-o) $(gen-f90)
	$(fc) -I$(lb) -J$(lb) $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=))



$(lb)/%.o:%.f90
	gfortran -I$(lb) -J$(lb) -c $< -o $@ -Wall

$(b)/%.o:%.f90
	$(compiler) -c $^ -o $@


$(lib-a):$(wildcard ../src/*.f90)
	(cd .. && make)

$(lb)/compressed_tensors.o: $(lb)/compressed_utils.o $(print-o) $(lb)/detrace_apple.o $(lb)/calc_derivs.o $(lb)/compressed_arrays.o
$(lb)/calc_derivs.o:$(lb)/sf_disp_tangtoe.o
$(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o: $(lb)/data_types.o

$(lb)/compressed_utils.o: $(print-o)

clean:
	rm *.x build/*

comp-clean:
	rm *.x build/comp*

