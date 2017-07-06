# Install SCME in QUIP

#? Alter this line to the QUIP root dir:
export QUIP_ROOT=~/gits/QUIP/

#1 Alter sourceMe_env.sh and link to QUIP/
ln -sf $PWD/sourceMe_quipenv.sh $QUIP_ROOT/sourceMe_quipenv.sh

#2 Create ThirdParty dir. Link/copy things into it. 
mkdir -f $QUIP_ROOT/src/ThirdParty/
ln -sf $PWD/cube_tools.f95 $QUIP_ROOT/src/ThirdParty/cube_tools.f95
ln -sf $PWD/Makefile.ThirdParty $QUIP_ROOT/src/ThirdParty/Makefile
ln -sf $PWD/../../scme $QUIP_ROOT/src/ThirdParty/scme


#3 Link IPModel_SCME.f95 in QUIP/src/Potenitals/ (substituting the present one)
ln -sf $PWD/IPModel_SCME.f95 $QUIP_ROOT/src/Potentials/IPModel_SCME.f95 

#4 Go to $QUIP_ROOT and source the environment, and install
(cd $QUIP_ROOT && . sourceMe_quipenv.sh)
#Only first time, to create Make.inc: (cd $QUIP_ROOT && make config) #only needed once 
(cd $QUIP_ROOT && make install-quippy)

