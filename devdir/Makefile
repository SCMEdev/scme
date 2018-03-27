#####################################################################
ext=.f90
bd := build
srcdirs=. ../src/

vpath %$(ext) $(srcdirs) 

sourcenames:= $(notdir $(wildcard $(srcdirs:=/*$(ext))))

objects:= $(addprefix $(bd)/, $(sourcenames:%$(ext)=%.o))


#########################################################################
#User Settings
fc := gfortran
opti:= 
#-Ofast -flto -march=native 
# -ftree-vectorize -ftree-loop-if-convert -ftree-loop-distribution -finline-functions 
err := -Wall -fmax-errors=3 -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow


base:= -I$(bd) -J$(bd) -ffree-line-length-0 -cpp

flags:= $(base) $(err) $(opti) 
# Must place main object first for generic program to call the correct module. 
generic = $^ -o $@ -cpp -DMODU=$(notdir $(<:.o=)) generic_program.f08

comp = $(fc) $(flags) $(generic)

vpath %.f90 ../src

# -I = where to search! 
# -J = where to put!    

all:
	make -j4 run.x
	make -j4 tests.x

$(bd)/%.o:%.f90
	@mkdir -p $(@D)
	$(fc) $(flags) -c $< -o $@

#apple detracer cannot run with bounds checking, so skip all the error flags. 
$(bd)/detrace_apple.o:detrace_apple.f90
	@mkdir -p $(@D)
	$(fc) $(base) $(opti) -c $< -o $@



#all:
#	make -j4 run.x
#	make -j4 tests.x

old: 
	make -j4 pold.x detraceo.x detraceh.x detraceg.x qpole.x opole.x

new: 
	make -j4 apple.x utils.x convert.x fool.x tests.x run.x 

#Executables
detraceo.x: $(bd)/detrace_mod.o $(objects)
	$(comp) -DSUB=io_octa


detraceh.x: $(bd)/detrace_mod.o $(objects)
	$(comp) -DSUB=io_hexa


detraceg.x: $(bd)/detrace_mod.o $(objects)
	$(comp) -DSUB=io_general


pold.x: $(bd)/poldip.o  $(objects)
	$(comp)

qpole.x: $(bd)/qpole_dev.o $(objects) 
	$(comp)

opole.x: $(bd)/opole_dev.o $(objects)
	$(comp)


### COMPRESSED

apple.x: $(bd)/detrace_apple.o 
	$(comp) 

# using 'objects':
utils.x: $(bd)/compressed_utils.o  $(objects)
	$(comp)
convert.x: $(bd)/compressed_tensors.o $(bd)/compressed_utils.o $(objects)
	$(comp)
fool.x: $(bd)/compressed_foolin.o $(objects) 
	$(comp)

tests.x: $(bd)/compressed_tests.o $(objects)
	$(comp)

run.x: $(bd)/compressed_run.o $(objects)
	$(comp)

clean:
	rm -rf $(bd) *.x


#####################################################################################


include prerequisites.makefile



