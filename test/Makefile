scme_obj = ../build
scme_mod = ../build
#test_build = build

FC = gfortran
FFLAGS = -O3 -cpp -pg -I$(scme_mod)


run.x: run_tests.o mifu_asserts.o test_scme.o
	$(FC) $(FFLAGS) -o run.x run_tests.o mifu_asserts.o test_scme.o  $(scme_obj)/libscme.a

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
