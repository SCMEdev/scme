from quippy import *
from ase import atoms
from ase.io import read,write
from ase.visualize import view
import os,sys
import numpy as np

# the quippy+ASE way of calculating energies

os.chdir(".")
hexamers = read("6mer.xyz")
#,index="1")
#pot=Potential("Potential xml_label=scme_plus_gap_h2o_2b_ccsdt_scme",param_filename="quip_params.xml")

#pot2b=Potential("Potential xml_label=scme_plus_gap_h2o_2b_ccsdt_scme",param_filename="gp.xml")
#pot2b=Potential("Potential xml_label=gap_h2o_2b_ccsdt_scme",param_filename="gp.xml")
#pot3b=Potential("Potential xml_label=GAP_2015_6_8_60_19_2_7_255",param_filename="gp_3b.xml")
#pot = Potential("Sum",pot1=pot2b,pot2=pot3b)

#pspot=Potential("IP PartridgeSchwenke",param_filename="ps.xml")



kcal=23.060541945
#ff=open('e.dat','w')
e2b = []
#for i in range(len(hexamers)):

## Finite Difference
print ""
print ""
print "Calculating Finite Difference Forces:"

pot=Potential("IP SCME force_using_fd=True label=version_20160315",param_filename="dummy.xml")

hexa = hexamers
hexa.set_cell([50.0,50.0,50.0])
hexa.set_pbc(True)

hexa.set_calculator(pot)
e = hexa.get_potential_energy()
f = hexa.get_forces()

s = hexa.get_chemical_symbols()
for i in range(len(s)/3):
	print "H2O Nr.: ", i+1
	print s[3*i +0],f[3*i+0,:]
	print s[3*i +1],f[3*i+1,:]
	print s[3*i +2],f[3*i+2,:]

print e

## Analytical forces!
print ""
print ""
print "... and now analytical forces:"

pot2=Potential("IP SCME label=version_20160315",param_filename="dummy.xml")

hexa2 = hexa

hexa2.set_cell([50.0,50.0,50.0])
hexa2.set_pbc(True)

hexa2.set_calculator(pot2)
e2 = hexa2.get_potential_energy()
f2 = hexa2.get_forces()

s2 = hexa2.get_chemical_symbols()
for i in range(len(s)/3):
	print "H2O Nr.: ", i+1
	print s2[3*i +0],f2[3*i+0,:]
	print s2[3*i +1],f2[3*i+1,:]
	print s2[3*i +2],f2[3*i+2,:]

print e2


