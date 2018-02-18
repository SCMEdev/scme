#settings

bd=build
srcdirs=src 

fflags:= $(opti) -pg -I$(bd) -J$(bd) -cpp -D'DEBUG_PRINTING' -ffree-line-length-0 -Wall
fc:=gfortran


target=$(bd)/libscme.a


#####################################################################
vpath %.f90 $(srcdirs) 

sourcenames:= $(notdir $(wildcard $(srcdirs:=/*.f90)))

objects:= $(addprefix $(bd)/, $(sourcenames:%.f90=%.o))

all:
	make -j4 $(target)
	

$(target): $(objects)
	ar rcs $@ $^

$(bd)/%.o: %.f90
	@mkdir -p $(bd)
	$(fc) $(fflags) -o $@ -c $<

clean:
	@rm -rf $(bd)

include prerequisites.makefile




