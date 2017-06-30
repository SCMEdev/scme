#Definitions
scme_obj = ../build
scme_mod = ../build
scme_src = ../src

#scme_sources = $(shell ls ../src/*.f90)
#scme_objects = $(patsubst  ../src/%, ../build/%, $(scme_sources:.f90=.o))

b=build


FC = gfortran
FFLAGS = -O0 -cpp -pg -I$(scme_mod) -J$b
#-fopenmp
## flags for home:
#-Ofast -cpp -pg -I$(scme_mod) -fopenmp
#-fopenmp 

#default
all:
	make scme
	make run.x
	make 2run.x

# commands:
klean_all__SCME_AND_TEST:
	(cd ../ && make clean) 
	make clean

diff_reference__RUN:
	./2run.x |diff - reference.dat

renew_reference_and_print__RUN:
	./2run.x |tee  reference.dat

just__RUN:
	./2run.x

clean:
	rm -f build/* *.x

clean-scme:
	rm -f ../build/* 


#executables:
common_obj=$(addprefix $b/, mifu_asserts.o test_scme.o )

run.x: $b/run_tests.o $(common_obj) $(scme_obj)/libscme.a


2run.x: $b/run_tests2.o $(common_obj) $(scme_obj)/libscme.a


#rules:
%.x:
	$(FC) $(FFLAGS) -o $@ $^ 

$b/%.o:%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$b:
	mkdir $@

scme:
	(cd ../ && make) 

#object dependencies:
$b/run_tests2.o\
$b/run_tests.o:\
$(common_obj)

$b/test_scme.o:\
$b/mifu_asserts.o



