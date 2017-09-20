
#give sourcefile on command line
ifndef f
abort:
	@echo "Call with 'make f=program_file.f90' "
endif

#scme paths/files
scme_root:=/home/j/gits/scme
scme_src:=$(scme_root)/src
scme_build:=$(scme_root)/build
scme_lib:=$(scme_build)/libscme.a

FC:=gfortran

local_src:=$f
local_obj:=$(local_src:.f90=.o)

prog:=$(local_src:.f90=.x)


#order matters here:
$(prog): $(local_obj) $(scme_lib) 
	$(FC)  $^ -o $@

%.o: %.f90 $(scme_lib)
	$(FC) -c $< -o $@ -I$(scme_build) 


$(scme_lib): $(wildcard $(scme_src)/*.f90)
	(cd $(scme_root) && make)


clean:
	rm -f *.o *.x *.mod

