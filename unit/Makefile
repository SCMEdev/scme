#Definitions
scme_root := ..
scme_build := $(scme_root)/build
scme_src := $(scme_root)/src

s := src
b := build


FC = gfortran
FFLAGS = -O0 -cpp -I$(scme_build) -J$b $(comline)
#-pg
#-fopenmp
## flags for home:
#-Ofast -cpp -pg -I$(scme_build) -fopenmp
#-fopenmp 

#default
all: $b
	make scme 
	make run.x
	make 2run.x

# commands:
klean_all__SCME_AND_TEST:
	(cd $(scme_root) && make clean) 
	make clean

diff_reference__RUN:
	./2run.x |diff - reference.dat

renew_reference_and_print__RUN:
	./2run.x |tee  reference.dat

just__RUN:
	./2run.x

clean:
	rm -rf build *.x $(scme_build)/* 


#executables:
common_obj=$(addprefix $b/, mifu_asserts.o test_scme.o ) 

run.x: $b/run_tests.o $(common_obj) $(scme_build)/libscme.a


2run.x: $b/run_tests2.o $(common_obj) $(scme_build)/libscme.a


#rules:
%.x:
	$(FC) $(FFLAGS) -o $@ $^ 

$b/%.o:$s/%.f90 
	$(FC) $(FFLAGS) -c $< -o $@

$b:
	mkdir $@

scme:
	(cd $(scme_root) && make) 

#object dependencies:
$b/run_tests2.o\
$b/run_tests.o:\
$(common_obj)

$b/test_scme.o:\
$b/mifu_asserts.o



