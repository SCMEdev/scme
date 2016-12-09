scme_obj = ../obj
scme_mod = ../mod
scme_src = ../src
#test_build = build

FC = gfortran
FFLAGS = -O2 -cpp -pg -I$(scme_mod) 
#-fopenmp 


run.x: run_tests.o mifu_asserts.o test_scme.o $(scme_obj)/libscme.a
	$(FC) $(FFLAGS) -o run.x $^


## The 9 and 15 H2O tests:
2run.x: run_tests2.o mifu_asserts.o test_scme.o $(scme_obj)/libscme.a
	$(FC) $(FFLAGS) -o $@ $^
	
run_tests2.o: run_tests2.f90 mifu_asserts.mod test_scme.mod
	$(FC) $(FFLAGS) -c run_tests2.f90
	
scme:
	(cd ../ && make) 

kscme:
	(cd ../ && make clean) 

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
