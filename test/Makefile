scme_obj = ../obj
scme_mod = ../mod
scme_src = ../src
#test_build = build

FC = gfortran
FFLAGS = -O0 -cpp -pg -I$(scme_mod) 
#-fopenmp
## flags for home:
#-Ofast -cpp -pg -I$(scme_mod) -fopenmp
#-fopenmp 

all:
	make scme
	make run.x
	make 2run.x

run.x: run_tests.o mifu_asserts.o test_scme.o $(scme_obj)/libscme.a
	$(FC) $(FFLAGS) -o run.x $^


## The 9 and 15 and 32 and 64 H2O tests:
2run.x: run_tests2.o mifu_asserts.o test_scme.o $(scme_obj)/libscme.a
	$(FC) $(FFLAGS) -o $@ $^

run_tests2.o: run_tests2.f90 mifu_asserts.mod test_scme.mod
	$(FC) $(FFLAGS) -c run_tests2.f90

scme:
	(cd ../ && make) 

klean_all__SCME_AND_TEST:
	(cd ../ && make clean) 
	make clean

diff_reference__RUN:
	./2run.x |diff - reference.dat

renew_reference_and_print__RUN:
	./2run.x |tee  reference.dat


#$(scme_obj)/libscme.a: $(wildcard ../$(scme_src)/*.f90)
#	cd ../ && make 
#run_tests.o mifu_asserts.o test_scme.o  $(scme_obj)/libscme.a



run_tests.o: run_tests.f90 mifu_asserts.mod test_scme.mod
	$(FC) $(FFLAGS) -c run_tests.f90

mifu_asserts.mod: mifu_asserts.o

mifu_asserts.o: mifu_asserts.f90
	$(FC) $(FFLAGS) -c mifu_asserts.f90

test_scme.mod: test_scme.o

test_scme.o: test_scme.f90 mifu_asserts.mod
	$(FC) $(FFLAGS) -c test_scme.f90

clean:
	rm -f *.o *.mod *.x *~
