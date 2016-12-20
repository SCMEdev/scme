dip: 
	gfortran ps_dip.f90 test_ps_dip.f90 -o pdip && ./pdip

pot: constants.o
	gfortran ps_pot.f90 test_ps_pot.f90 constants.o -o ppot 
	./ppot

potr: 
	gfortran mod_ps_potr.f90 -o ppotr && ./ppotr

cpp_dip: ps.o constants.h ps.h
	gfortran test_ps_cpp.f90 ps.o -o cppdip 
	./cppdip

lodi_dip:
	gfortran lodi_dip.f90 test_lodi_dip.f90 -o plodi 
	./plodi

dms: constants.o
	gfortran ps_dms.f90 test_ps_dms.f90 constants.o -o pdms 
	./pdms

cpp_pot: ps.o constants.h ps.h
	gfortran test_cpp_pot.f90 ps.o -o ppot_cpp 
	./ppot_cpp


constants.o: constants.f90
	gfortran -c ../src/data_types.f90 constants.f90  

	
ps.o: ps.cpp ps.h
	g++ -c ps.cpp

clean: 
	rm ppot pdip *.mod

