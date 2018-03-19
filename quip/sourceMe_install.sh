
# ** SETTINGS **

# Alter these lines to what you use when installing QUIP
export QUIP_ROOT=~/gits/QUIP/
export QUIP_ARCH=linux_x86_64_gfortran_openmp


echo "
export QUIP_ARCH=$QUIP_ARCH
export QUIP_ROOT=$QUIP_ROOT
# . activate my_python2.7_virtualenv
" > $QUIP_ROOT/sourceMe_environment.sh




#########################################################################
#########################################################################
#########################################################################

# ** DONT TOUCH (UNLESS YOU WANT TO) **
# Create/keep ThirdParty dir. Link things into it. 
mkdir -pv $QUIP_ROOT/src/ThirdParty/
ln -sfv $PWD/cube_tools.f95 		$QUIP_ROOT/src/ThirdParty/cube_tools.f95
ln -sfv $PWD/Makefile.ThirdParty $QUIP_ROOT/src/ThirdParty/Makefile
mv -v $QUIP_ROOT/src/ThirdParty/scme  $QUIP_ROOT/src/ThirdParty/scme_BACKUP_$(date +"%m-%d-%Y_%H:%M:%S")
ln -sfv $PWD/../../scme 			$QUIP_ROOT/src/ThirdParty/scme

# Link IPModel_SCME.f95 in QUIP/src/Potenitals/ 
ln -sfv $PWD/IPModel_SCME.f95 $QUIP_ROOT/src/Potentials/IPModel_SCME.f95 

echo ""
echo "# Done. 
# Go to $QUIP_ROOT and do: 

. sourceMe_environment.sh
make config
make install-quippy
"


