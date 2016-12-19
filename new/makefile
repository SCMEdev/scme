dip: 
	gfortran ps_dip.f90 test_ps_dip.f90 -o pdip && ./pdip

pot: 
	gfortran ps_pot.f90 test_ps_pot.f90 -o ppot && ./ppot

potr: 
	gfortran mod_ps_potr.f90 -o ppotr && ./ppotr

cpp_dip_test:
	g++ -c ps.cpp && gfortran -c test_ps_cpp.f90 && gfortran ps.o test_ps_cpp.o -o cppdip && ./cppdip

lodi_dip_test:
	gfortran lodi_dip.f90 test_lodi_dip.f90 -o plodi && ./plodi

clean: 
	rm ppot pdip *.mod

