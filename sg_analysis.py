#!/usr/bin/env python3

# just a simple script to retrieve some information on the geometry and symmetry of a given POSCAR, CONTCAR, or CIF-file

from pymatgen.core import Composition,Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys

name = sys.argv[1]
struc = Structure.from_file(name)


composition = Composition(struc.composition)

formula = composition.reduced_formula
#print(composition.get_reduced_composition_and_factor()[1])

# analyze space group
space_group = SpacegroupAnalyzer(struc)

HM_number = space_group.get_space_group_number()
HM_desc = space_group.get_space_group_symbol()
Z = composition.get_reduced_composition_and_factor()[1]
V = struc.volume

# calculate cubic optimum
a_opt = V**(1/3)
latt_opt = 0.0

print("The lattice parameters a, b, and c in Angstrom are:")
for x in struc.lattice.abc:
    print("{:8.3f}".format(x))
    latt_opt += (x-a_opt)**2
    
latt_opt = (latt_opt)**(1/2)


print('The space group of this structure is '+HM_desc+', the corresponding HM number is '+str(HM_number)+".\nThe lattice contains "+str(Z)+" formula unit(s).")
print('The structure has a volume of '+"{:10.3f}".format(V)+".")
print('Divergence from cube:',"{:6.3f}".format(latt_opt/a_opt))

# tra = space_group.get_conventional_to_primitive_transformation_matrix
# print(tra)


