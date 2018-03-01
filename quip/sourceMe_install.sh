
# ** SETTINGS **
# Alter this line to the QUIP root dir that you have cloned:
export QUIP_ROOT=~/gits/QUIP/

# Alter sourceMe_env.sh and link to QUIP/, since it will be convenient to have it there when recompiling quippy. 
ln -sf $PWD/sourceMe_quipenv.sh $QUIP_ROOT/sourceMe_quipenv.sh


# ** LINK SOURCE FILES **
# Create/keep ThirdParty dir. Link things into it. 
mkdir -p $QUIP_ROOT/src/ThirdParty/
ln -sf $PWD/cube_tools.f95 		$QUIP_ROOT/src/ThirdParty/cube_tools.f95
ln -sf $PWD/Makefile.ThirdParty $QUIP_ROOT/src/ThirdParty/Makefile
mv  $QUIP_ROOT/src/ThirdParty/scme  $QUIP_ROOT/src/ThirdParty/scme_BACKUP_$(date +"%m-%d-%Y_%H:%M:%S")
ln -sf $PWD/../../scme 			$QUIP_ROOT/src/ThirdParty/scme

# Link IPModel_SCME.f95 in QUIP/src/Potenitals/ 
ln -sf $PWD/IPModel_SCME.f95 $QUIP_ROOT/src/Potentials/IPModel_SCME.f95 

echo ""
echo "Installed SCME into QUIP!

go to QUIP/ and do 
. sourceMe_quipenv.sh
make config
make install-quippy

or from here

. sourceMe_quipenv.sh
(cd $QUIP_ROOT && make config) #only needed the first time! 
(cd $QUIP_ROOT && make install-quippy)
"


