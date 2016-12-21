my_older_ps_dip: 
	gfortran ps_dip.f90 test_ps_dip.f90 -o pdip && ./pdip

lodi_dip:
	gfortran lodi_dip.f90 test_lodi_dip.f90 -o plodi 
	./plodi

potr: 
	gfortran mod_ps_potr.f90 -o ppotr && ./ppotr

#///////////////////////////////////////////////////////////////////////
pot: constants.o
	gfortran ps_pot.f90 test_ps_pot.f90 constants.o -o ppot 
	./ppot

dms: constants.o
	gfortran ps_dms.f90 test_ps_dms.f90 constants.o -o pdms 
	./pdms

cpp_dip: ps.o constants.h ps.h
	gfortran test_cpp_dms.f90 ps.o -o cppdip 
	./cppdip

cpp_dip_src: ../obj/ps.o ../src/constants.h ../src/ps.h
	gfortran test_cpp_dms.f90 ../obj/ps.o -o cppdip 
	./cppdip


cpp_pot: ps.o constants.h ps.h
	gfortran test_cpp_pot.f90 ps.o -o ppot_cpp 
	./ppot_cpp
#///////////////////////////////////////////////////////////////////////
pconst: test_constants.f90 constants.o 
	gfortran test_constants.f90 constants.o -o pconst
	./pconst
#///////////////////////////////////////////////////////////////////////
constants.o: constants.f90
	gfortran -c ../src/data_types.f90 constants.f90  

	
ps.o: ps.cpp ps.h
	g++ -c ps.cpp

constants.h: ../src/constants.h
	cp $< . 

ps.h: ../src/ps.h
	cp $< .

#///////////////////////////////////////////////////////////////////////
clean: 
	rm ppot pdip *.mod

