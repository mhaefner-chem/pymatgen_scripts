#!/usr/bin/env python3

# just a simple script to retrieve some information on the geometry and symmetry of a given POSCAR, CONTCAR, or CIF-file

from pymatgen.core import Composition,Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import sys

def pos_extract(str):
    pos = float(str.split("(")[0])
    return pos


sys.stdout.flush()

max_p = 0.25
max_sc = 1

n = len(sys.argv)
for i in range(0,n):
    # specify the cell
    if sys.argv[i] == "-c":
        print("Reading in the file "+sys.argv[i+1]+".")
        struc = Structure.from_file(sys.argv[i+1])
        name = sys.argv[i+1]
        
    # dopant element + oxi state
    if sys.argv[i] == "-d":
        dope_el = sys.argv[i+1]
        dope_oxi = int(sys.argv[i+2])
    
    # doped element
    if sys.argv[i] == "-r":
        repl = sys.argv[i+1]
        
    # maximum doping percentage
    if sys.argv[i] == "-p":
        max_p = float(sys.argv[i+1])
        
    # maximum supercell size
    if sys.argv[i] == "-s":
        max_sc = int(sys.argv[i+1])
        
    
        
        
        


elements = ["H","He","Li","Be","B","C","O","N","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br"]


type = []
site = []
oxi_state = {}
occupancy = {}
pos = {}
pos_single = [0,0,0]
multi = {}
site_type = {}
split_lines = []
sites = []


stopper = ["loop_","End","\#"]
switch = 0


with open(name,'r') as file:
    lines = file.readlines()


    for i in range(0,len(lines)):
        
        split_lines.append(lines[i].split())
        
        if "_cell_formula_units_Z" in lines[i]:
            Z = lines[i].split()[1]        
        if "_atom_type" in lines[i]:
            type.append(lines[i])
            switch = 1
            if "_atom_type_symbol" in lines[i]:
                type_symbol = len(type)-1
            if "_atom_type_oxidation_number" in lines[i]:
                type_oxi_num = len(type)-1
        
        if "_atom_site" in lines[i]:
            site.append(lines[i])
            switch = 2
            if "_atom_site_occupancy" in lines[i]:
                site_occup = len(site)-1
            if "_atom_site_type_symbol" in lines[i]:
                site_type_symbol = len(site)-1
            if "_atom_site_label" in lines[i]:
                site_label = len(site)-1
            if "_atom_site_symmetry_multiplicity" in lines[i]:
                site_multi = len(site)-1
            if "_atom_site_Wyckoff_symbol" in lines[i]:
                site_Wyckoff = len(site)-1
            if "_atom_site_fract_x" in lines[i]:
                site_x = len(site)-1
            if "_atom_site_fract_y" in lines[i]:
                site_y = len(site)-1
            if "_atom_site_fract_z" in lines[i]:
                site_z = len(site)-1
                
                
        if not "_" in lines[i]:
            if switch == 1:
                for j in range(i,len(lines)):
                    if any(item in lines[j] for item in stopper):
                        switch = 0
                        break
                    print(lines[j])
                    oxi_state[lines[j].split()[type_symbol]] = int(lines[j].split()[type_oxi_num])
                    
            if switch == 2:
                if any(item in lines[i] for item in stopper):
                    switch = 0
                    break

                key = lines[i].split()[site_label]
                occupancy[key] = float(lines[i].split()[site_occup])
                multi[key] = int(lines[i].split()[site_multi])
                site_type[key] = lines[i].split()[site_type_symbol]
                sites.append(key)
                

                pos_single[0] = pos_extract(lines[i].split()[site_x])
                pos_single[1] = pos_extract(lines[i].split()[site_y])
                pos_single[2] = pos_extract(lines[i].split()[site_z])
                
                pos[key] = [pos_single[0],pos_single[1],pos_single[2]]
                print(key,pos_single)
                

repl_unit = 0
for site in sites:
    if repl in site:
        repl_unit += multi[site]*occupancy[site]

repl_tot = round(max_sc*repl_unit)

for i in range(1,repl_tot):
    if i/repl_tot > max_p:
        break
    print("{:6.3f}".format(i/repl_tot*100.0),"%")
    
    
    



                # split_lines[i][site_type_symbol] = "XYZ"
                    



    # for i in range(j,len(lines)):
        
    

            # for j in range(i,len(lines)):
                
    #     sys.exit()
            
        
    #     if 'E_0' in line:
    #         ref_e[ref_name[reference]] = line.split()[1]
    #     if 'T/K' in line:
    #         switch = 1
    #         continue
    #     if switch == 1:
    #         thermo.append(line.split()[1])
    # ref_thermo[ref_name[reference]] = thermo

sys.exit()


composition = Composition(struc.composition)

formula = composition.reduced_formula
#print(composition.get_reduced_composition_and_factor()[1])

# analyze space group
space_group = SpacegroupAnalyzer(struc)

HM_number = space_group.get_space_group_number()
HM_desc = space_group.get_space_group_symbol()
Z = composition.get_reduced_composition_and_factor()[1]
V = struc.volume

# for site in struc.sites:
#     print(site.species)
    
# transformation = dope("Al3+",alio_tol=2,codopant=True,allowed_doping_species=["Al3+","Na1+","Cl1-"])
# doped_struc = transformation.apply_transformation(structure=struc) 
# doped_struc.to("test.cif",fmt="cif")


# a = struc.apply_transformation("Al3+",alio_tol=1,codopant=True,allowed_doping_species=["Na1+"])


print("The lattice parameters a, b, and c in Angstrom are:")
for x in struc.lattice.abc:
    print("{:8.3f}".format(x))

print('The space group of this structure is '+HM_desc+', the corresponding HM number is '+str(HM_number)+".\nThe lattice contains "+str(Z)+" formula unit(s).")
print('The structure has a volume of '+str(V)+".")
