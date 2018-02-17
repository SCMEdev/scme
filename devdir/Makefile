### DIRECTORIES
lb := lbuild
dodir=@mkdir -p $(@D)


# Fast commands
all:
	make -j4 run.x
	make -j4 tests.x

old: 
	make -j4 pold.x detraceo.x detraceh.x detraceg.x qpole.x opole.x

new: 
	make -j4 apple.x utils.x convert.x fool.x tests.x run.x 



### COMPILATION 
fc := gfortran
opti:= 
#-Ofast -flto -march=native 
# -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -finline-functions 
err := -Wall -fmax-errors=3 -fcheck=all

flags:= -I$(lb) -J$(lb) $(err) $(opti) -ffree-line-length-0

# Must place main object first for generic program to call the correct module. 
generic = $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) generic_program.f90

comp2 = $(fc) $(flags) $(generic)

vpath %.f90 ../src

# -I = where to search! 
# -J = where to put!    


$(lb)/%.o:%.f90
	$(dodir)
	$(fc) $(flags) -c $< -o $@


#$(lb)/detrace_apple.o:detrace_apple.f90
#	$(dodir)
#	$(fc) -I$(lb) -J$(lb) -c $< -o $@


scme_objs0:=\
$(lb)/localAxes_mod.o\
$(lb)/opole.o\
$(lb)/qpole.o\
$(lb)/polariz_parameters.o\
$(lb)/localAxes_mod.o\
$(lb)/multipole_parameters.o\
$(lb)/sf_disp_tangtoe.o\


scme_objs:=\
$(scme_objs0)\
$(lb)/printer_mod.o\
$(lb)/data_types.o\

# SCME DEPENDENCIES

$(lb)/qpole_dev.o\
$(lb)/opole_dev.o\
$(lb)/detrace_mod.o:\
$(scme_objs)

$(scme_objs0):$(lb)/printer_mod.o
$(lb)/printer_mod.o:$(lb)/data_types.o
$(lb)/calc_derivs.o:$(lb)/sf_disp_tangtoe.o


$(lb)/opole.o:$(lb)/printer_mod.o 


detraceo.x: $(lb)/detrace_mod.o $(scme_objs)
	$(comp2) -DSUB=io_octa


detraceh.x: $(lb)/detrace_mod.o $(scme_objs)
	$(comp2) -DSUB=io_hexa


detraceg.x: $(lb)/detrace_mod.o $(scme_objs)
	$(comp2) -DSUB=io_general


pold.x: $(lb)/poldip.o  $(scme_objs) 
	$(comp2)

qpole.x: $(lb)/qpole_dev.o $(scme_objs) 
	$(comp2)

opole.x: $(lb)/opole_dev.o $(scme_objs) 
	$(comp2)







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

test_objects := $(objects) $(lb)/compressed_foolin.o $(lb)/compressed_tensors.o $(lb)/compressed_utils.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o $(lb)/sf_disp_tangtoe.o $(lb)/data_types.o $(lb)/polariz_parameters.o
# using 'test_objects':
tests.x: $(lb)/compressed_tests.o $(test_objects)
	$(comp2)

run.x: $(lb)/compressed_run.o $(lb)/compressed_tests.o $(test_objects)
	$(comp2)


### DEPENDENCIES

$(lb)/compressed_utils.o: $(lb)/printer_mod.o $(lb)/compressed_arrays.o

$(lb)/compressed_tensors.o: $(lb)/compressed_utils.o 

$(lb)/compressed_foolin.o: $(lb)/compressed_tensors.o 

$(lb)/compressed_tests.o: $(lb)/compressed_tensors.o $(lb)/detrace_apple.o $(lb)/calc_derivs.o $(lb)/polariz_parameters.o

$(lb)/compressed_run.o: $(lb)/compressed_tests.o 




### CLEAN
clean:
	rm -rf *.x *.mod *.o lbuild/ obuild/

comp-clean:
	rm *.x lbuild/comp*

