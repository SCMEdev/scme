# To install SCME into QUIP and install quippy into your desired python 2.7 (virtual) environment:
# 0. Place the scme/ in the same directory as QUIP (NOT inside QUIP/src/Thirdparty)
# 1. Alter the QUIP_ROOT path in both "sourceMe_" files.
# 2. Alter the ". activate ..." in 'sourceMe_quipenv.sh' to point to your virtualenv, or remove it if you want default installation place.
# 3. Do ". sourceMe_install.sh" and follow the instrucitons on how to compile. 
# 4. Go to scme/test/quippy and run the python script to see if it works. 
# 
# 
# These are some basic example commands to use the scme potential using quippy in python 2.7: 
# 


# NB. quippy only works in python 2

import numpy as np
from quippy import Potential
from ase import atoms
from ase.io import read,write
from ase.visualize import view
import os,sys


scmepot = Potential("IP SCME label=version_20170802",param_string="",
          use_repulsion=False,use_super_repulsion=True,full_interaction_order=True,use_variable_quadrupole=True,
          use_variable_octapole=False, use_compressed=True,damping_parameter=1.2,kernel_choice=1)


structure = read("temp.xyz") #change to existing file

structure.set_cell([50.0,50.0,50.0])
structure.set_pbc(True)

structure.set_calculator(scmepot)

print structure.get_potential_energy()

# structure.
# structure.
# structure.
