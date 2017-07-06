## Variables
prog:=prog.x
b:=../../build
s:=../../src
F:=gfortran
flags= -I$b -J$b -cpp -D$@
common:= $b/data_types.o $b/localAxes_mod.o

## Compilation
laxes:program.f90
	$F $(flags)  $(common) $b/printer_mod.o $b/ps_dms.o $< -o $(prog)

gdamp:program.f90
	$F $(flags) $(common) $< -o $(prog)


## Commands
run:
	./$(prog)

clean: 
	rm *.mod *.o a.out



#check:
#ifndef dothis
#dothis:=gdamp
#else ifeq ($(dothis),)
#dothis:=gdamp
#else ifeq ($(dothis),"")
#dothis:=gdamp
#endif

	
#dothis?=gdamp
#checkenv:
#ifndef DOTHIS
#DOTHIS:=gdamp
#endif 

#	$(error "DOTHIS not defined")
#DOTHIS:=gdamp
