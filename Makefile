fflags:= $(opti) -pg -I$(bd) -J$(bd) -cpp -D'DEBUG_PRINTING' -ffree-line-length-0 -Wall
fc:=gfortran

bd:=build
srcdirs:=src devdir
ext:=f90


#Settings above. Furst run './prescan src devdir' to create the dependency file 'prerequisites.makefile'
#####################################################################
#####################################################################
vpath %.f90 $(srcdirs) 

srcs:= $(notdir $(wildcard $(srcdirs:=/*.$(ext))))

objects:= $(addprefix $(bd)/, $(srcs:%.$(ext)=%.o))

#####################################################################

all:
	make -j4 target


target: $(objects)
	@echo "objects are in $(bd) !!!"

$(bd)/%.o: %.f90
	@mkdir -p $(bd)
	$(fc) $(fflags) -o $@ -c $<

clean:
	rm -rf $(bd)

#####################################################################

include prerequisites.makefile




