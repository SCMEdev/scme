#####################################################################
bd=build
srcdirs=src 


vpath %.f90 $(srcdirs) 

sourcenames:= $(notdir $(wildcard $(srcdirs:=/*.f90)))

objects:= $(addprefix $(bd)/, $(sourcenames:%.f90=%.o))

target=$(bd)/libscme.a

#####################################################################

fflags:= $(opti) -pg -I$(bd) -J$(bd) -cpp -D'DEBUG_PRINTING' -ffree-line-length-0 -Wall
fc:=gfortran

all:
	make -j4 $(target)


$(target): $(objects)
	@echo "objects are in $(bd) !!!"

$(bd)/%.o: %.f90
	@mkdir -p $(bd)
	$(fc) $(fflags) -o $@ -c $<

clean:
	@rm -rf $(bd)

#####################################################################
include prerequisites.makefile




