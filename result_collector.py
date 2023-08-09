#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:40:47 2023

@author: bt308570
"""

# this program collects all results and compiles them into one file for each structure

import sys,os,importlib,glob
from pymatgen.core import Structure,Composition


sys.stdout.flush()
print("Collecting Results.")

# locate head directory
while True:
    if not os.path.isfile(".anchor"):
        os.chdir("..")
    else:
        print("Head directory found.")
        workdir = os.getcwd()
        break
    if os.getcwd() == "/":
        print("No anchor found.")
        sys.exit()
        
spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
directory_tools = importlib.util.module_from_spec(spec)
sys.modules["name"] = directory_tools
spec.loader.exec_module(directory_tools)        

 # move through the results
res_dir = workdir+"/RESULTS"
results_preliminary = [f for f in os.listdir(res_dir) if os.path.isdir(os.path.join(res_dir, f))]

results = []
for result in results_preliminary:
    if "ML" in result:
        continue
    else:
        results.append(result)
os.chdir(res_dir)


result_e = {}
result_thermo = {}
result_comp = {}

for result in results:
    thermo = []
    
    # create file with results, add name of the structure
    os.chdir(result)
    res_file = res_dir+"/"+result+'.txt'
    with open(res_file, 'w') as file:
        file.write(result+"\n")
        
    # if there's a file with energies, read the energy/unit
    if os.path.isfile("energy.txt"):
        with open('energy.txt','r') as file:
            for line in file:
                if 'Energy/unit' in line:
                    e_unit=(line.split()[1])
                    result_e[result] = e_unit
    
        # read the structure file and extract information about structure
        if os.path.isfile(workdir+"/OPT/"+result+"/CONTCAR"):
            structure = Structure.from_file(workdir+"/OPT/"+result+"/CONTCAR")
            atoms = structure.num_sites
            composition = Composition(structure.composition)
            result_comp[result] = composition
            frac = {}
            for element in composition.elements:
                frac[element] = composition.get_atomic_fraction(element)
        else:
            atoms = -1
                    
        # append structure composition, electronic energy, and atoms/cell    
        with open(res_file, 'a') as file:
            for i in frac:
                file.write("        "+i.symbol+"    "+"{:5.3f}".format(frac[i])+"\n")
            file.write("    E_0/eV "+"{:8.3f}".format(float(e_unit))+"\n")
            file.write("atoms/cell "+"{:8}".format(atoms)+"\n")
            
    # if there's a file with thermodynamic quantities, read the temperature and free energy/unit   
    temperatures = []
    free_energies_unit = []
    if os.path.isfile("thermo.txt"):
        
        with open('thermo.txt','r') as file:
            for line in file:
                if not 'unit' in line:
                    temperatures.append(line.split()[0])  
                    free_energies_unit.append(line.split()[2])  
    else:
        print("No thermodynamics for:",result)
        temperatures = [-100]
        free_energies_unit = [100]
                    
    # append the temperature and free energy/unit
    with open(res_file, 'a') as file:
        file.write("       T/K     E/eV\n")
        for i in range(0,len(temperatures)):
            t = "{:10.2f}".format(float(temperatures[i]))
            g_unit = "{:8.3f}".format(float(free_energies_unit[i]))
            file.write(t+' '+g_unit+"\n")
            thermo.append(g_unit)
    result_thermo[result] = thermo

        
    os.chdir(res_dir)
    
    
# find reference states in Settings

references = glob.glob(workdir+"/SETTINGS/*.vasp")
ref_struc = {}
ref_name = {}
ref_e = {}
ref_thermo = {}
ref_comp = {}
ref_elements = {}


for reference in references:
    ref_struc[reference] = Structure.from_file(reference)
    ref_name[reference] = directory_tools.directory_name_constructor(ref_struc[reference])
    with open(res_dir+'/'+ref_name[reference]+'.txt','r') as file:
        lines = file.readlines()
        switch = 0
        thermo = []
        for i in range(0,len(lines)):
            line = lines[i]
            if 'E_0' in line:
                ref_e[ref_name[reference]] = line.split()[1]
            if 'T/K' in line:
                switch = 1
                continue
            if switch == 1:
                thermo.append(line.split()[1])
        ref_thermo[ref_name[reference]] = thermo
    ref_comp[ref_name[reference]] = Composition(ref_struc[reference].composition)
    
    Z = ref_comp[ref_name[reference]].get_reduced_composition_and_factor()[1]
    num_atoms = ref_comp[ref_name[reference]].num_atoms
    ref_element = {}
    for element in ref_comp[ref_name[reference]]:
        ref_element[element.symbol] = ref_comp[ref_name[reference]].get_atomic_fraction(element)*num_atoms/Z
    ref_elements[ref_name[reference]] = ref_element

# plug everything together

for result in results:
    
    os.chdir(result)
    temperatures = []
    if os.path.isfile("thermo.txt"):
        
        with open('thermo.txt','r') as file:
            for line in file:
                if not 'unit' in line:
                    temperatures.append(line.split()[0])
    else:
        temperatures = [-100]
    os.chdir(res_dir)
    
    Z = result_comp[result].get_reduced_composition_and_factor()[1]
    num_atoms = result_comp[result].num_atoms
    
    res_elements = {}
    for element in result_comp[result]:
        res_elements[element.symbol] = result_comp[result].get_atomic_fraction(element)*num_atoms/Z
    
    share = {}
    
    delta_e = float(result_e[result])
    delta_thermo = []
    for i in range(0,len(temperatures)):
        # if 
        delta_thermo.append(float(result_thermo[result][i]))
    
    delta_file = res_dir+"/"+result+'_delta.txt'
    with open(delta_file, 'w') as file:
        file.write(result+"\n")
    
    frac_tot = {}
    
    for reference in references:
        share[ref_name[reference]] = 0
        frac = []
        for i in ref_elements[ref_name[reference]]:
            if i in res_elements:
                frac.append(res_elements[i]/ref_elements[ref_name[reference]][i])
            else:
                frac = [0]
                break
        
            
        #print(min(frac),ref_name[reference],result)
    
        delta_e = delta_e - (min(frac) * float(ref_e[ref_name[reference]]))
        for i in range(0,len(temperatures)):
            delta_thermo[i] = delta_thermo[i] - (min(frac) * float(ref_thermo[ref_name[reference]][i]))
            frac_tot[reference] = min(frac)
     
    tot = 0
    for fraction in frac_tot:
        tot += frac_tot[fraction]
        
    if tot == 0:
        print(result)
        tot = -1
            
    for fraction in frac_tot:
        with open(delta_file, 'a') as file:
            file.write("{:16}".format(ref_name[fraction])+" "+"{:4.0f}".format(frac_tot[fraction])+"{:7.2f}".format(frac_tot[fraction]/tot*100.0)+" %\n")
            
       
    for i in range(0,len(temperatures)):
        delta_thermo[i] = delta_thermo[i] + delta_e
        #print("{:10.2f}".format(float(temperatures[i])),"{:8.3f}".format(delta_thermo[i]))    
        
    with open(delta_file, 'a') as file:
        file.write("atoms/cell   "+"{:8.0f}".format(num_atoms)+"\n")
        file.write("Delta_E_0/unit "+"{:10.4f}".format(delta_e)+"\n")
        file.write("Delta_E_0/atom "+"{:10.4f}".format(delta_e/(num_atoms/Z))+"\n")

        file.write("       T/K       ΔE/u     ΔE/a\n")
        for i in range(0,len(temperatures)):
            t = "{:10.2f}".format(float(temperatures[i]))
            g_unit = "{:8.4f}".format(float(delta_thermo[i]))
            g_atom = "{:8.4f}".format(float(delta_thermo[i]/(num_atoms/Z)))
            file.write(t+'   '+g_unit+" "+g_atom+"\n")
    
            
#             #print(composition.get_atomic_fraction(element),element.symbol)
#             ele_ref_ratio = ref_comp[ref_name[reference]].get_atomic_fraction(element)
#             ele_res_ratio = result_comp[result].get_atomic_fraction(element)
            
#             Z = composition.get_reduced_composition_and_factor()[1]
#             print(result_comp[result].num_atoms)
#             if ele_res_ratio/ele_ref_ratio > ratio:
#                 ratio = ele_res_ratio/ele_ref_ratio
#         print(ratio,ref_comp[ref_name[reference]].reduced_formula,result_comp[result].reduced_formula)
