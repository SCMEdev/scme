dip: 
	gfortran ps_dip.f90 test_ps_dip.f90 -o pdip && ./pdip

pot: 
	gfortran ps_pot.f90 test_ps_pot.f90 -o ppot && ./ppot

potr: 
	gfortran mod_ps_potr.f90 -o ppotr && ./ppotr

cpp_dip_test:
	g++ -Ofast -c ps.cpp && gfortran -Ofast -c test_ps_cpp.f90 && gfortran -Ofast ps.o test_ps_cpp.o -o cppdip && ./cppdip

lodi_dip_test:
	gfortran lodi_dip.f90 test_lodi_dip.f90 -o plodi && ./plodi

dms: constants.o
	gfortran -Ofast ps_dms.f90 test_ps_dms.f90 constants.o -o pdms && ./pdms

constants.o:
	gfortran -c ../src/data_types.f90 constants.f90  

clean: 
	rm ppot pdip *.mod

