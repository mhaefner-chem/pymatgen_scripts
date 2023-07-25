#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:40:47 2023

@author: bt308570
"""

# this program collects all results and compiles them into one file for each structure

import sys,os,importlib,glob,warnings
from pymatgen.core import Structure,Composition
from pymatgen.io.vasp.outputs import Vasprun



def energy(outcar):
    last_match = None

    
    with open(outcar, 'r') as file:
        for line in file:
            if "free  e" in line:
                last_match = line
        
    last_match = last_match.split("=")[1]
    last_match = last_match.split("eV")[0]
    return float(last_match)

# beginning of main program
sys.stdout.flush()
warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
warnings.filterwarnings("ignore", message="your POTCAR is corrupted")
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

 # move into the preopt folder, gather result folders
preopt_dir = workdir+"/SUPERCELL/PREOPT"
os.chdir(preopt_dir)
dirs = [dir for dir in os.listdir(preopt_dir) if os.path.isdir(os.path.join(preopt_dir, dir))]

ref_energy = {}
ref_struc = {}
ref_comp = {}
references =  []
res_energy = {}
res_struc = {}
res_comp = {}
results = []

for dir in dirs:
    os.chdir(dir)
    if "base" in dir:
        
        print("Found a supercell structure",dir)
        
        configs = [config for config in os.listdir(os.path.join(preopt_dir,dir)) if os.path.isdir(os.path.join(preopt_dir,dir,config))]
        
        config_energy = {}

        
        for config in configs:
            os.chdir(config)
            if os.path.isfile("done"):
                struc = Structure.from_file("CONTCAR")
                composition = Composition(struc.composition)
                Z = composition.get_reduced_composition_and_factor()[1]        
                
                final_energy = energy("OUTCAR")
                config_energy[config] = final_energy/Z
                print(config,"{:8.3f}".format(final_energy/Z))
            else:
                config_energy[config] = -100.0
                if os.path.isfile("SCFCONV/done"):
                    xml = Vasprun("SCFCONV/vasprun.xml")
                    if xml.converged_electronic == False:
                        print("The SCF convergence failed!")
                else:
                    print(config,"Config not done!")
                
            os.chdir(os.path.join(preopt_dir,dir))
        res_energy[dir] = config_energy
        res_struc[dir] = struc
        res_comp[dir] = composition
        results.append(dir)

        
        
    else:
        print("Found a reference structure",dir)
        if os.path.isfile("done"):
            struc = Structure.from_file("CONTCAR")
            composition = Composition(struc.composition)
            Z = composition.get_reduced_composition_and_factor()[1]        
            
            final_energy = energy("OUTCAR")
            print("{:8.3f}".format(final_energy/Z))
        else:
             final_energy = -100.0
             
        
        ref_energy[dir] = final_energy/Z
        ref_struc[dir] = struc
        ref_comp[dir] = composition
        references.append(dir)

    os.chdir(preopt_dir)
    
    
sys.exit()    

# PROGRAM END REACHED



# evaluation with reference systems doesn't make sense yet

ref_elements = {}
for reference in references:
    
    
    Z = ref_comp[reference].get_reduced_composition_and_factor()[1]
    num_atoms = ref_comp[reference].num_atoms
    
    elements = {}
    for element in ref_comp[reference]:
        elements[element.symbol] = ref_comp[reference].get_atomic_fraction(element)*num_atoms/Z
    ref_elements[reference] = elements

    
for result in results:
    
    os.chdir(dir)
    res_elements = {}
    Z = res_comp[result].get_reduced_composition_and_factor()[1]
    num_atoms = res_comp[result].num_atoms
    

    for element in res_comp[result]:
        res_elements[element.symbol] = res_comp[result].get_atomic_fraction(element)*num_atoms/Z
        
    ref = Structure.from_file("ref.cif")
    ref_name = directory_tools.directory_name_constructor(ref)
    dope = Structure.from_file("dope.cif")
    dope_name = directory_tools.directory_name_constructor(dope)
    
    for reference in references:
        frac = {}
        
        if reference == ref_name or reference == dope_name:
            for i in ref_elements[reference]:
                if i in res_elements:
                    frac.append(res_elements[i]/ref_elements[reference][i])
                else:
                    frac = [0]
                    break
    os.chdir(preopt_dir)






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
            
        with open(delta_file, 'a') as file:
            file.write("{:16}".format(ref_name[reference])+" "+"{:4.0f}".format(min(frac))+"\n")
            
       
    for i in range(0,len(temperatures)):
        delta_thermo[i] = delta_thermo[i] + delta_e
        delta_thermo[i] = delta_thermo[i]/(num_atoms/Z)
        #print("{:10.2f}".format(float(temperatures[i])),"{:8.3f}".format(delta_thermo[i]))    
        
    with open(delta_file, 'a') as file:
        file.write("atoms/cell   "+"{:8.0f}".format(num_atoms)+"\n")
        file.write("Delta_E_0/atom "+"{:10.3f}".format(delta_e/(num_atoms/Z))+"\n")

        file.write("       T/K   Delta E/atom\n")
        for i in range(0,len(temperatures)):
            t = "{:10.2f}".format(float(temperatures[i]))
            g_unit = "{:8.3f}".format(float(delta_thermo[i]))
            file.write(t+'   '+g_unit+"\n")
    
            
#             #print(composition.get_atomic_fraction(element),element.symbol)
#             ele_ref_ratio = ref_comp[ref_name[reference]].get_atomic_fraction(element)
#             ele_res_ratio = result_comp[result].get_atomic_fraction(element)
            
#             Z = composition.get_reduced_composition_and_factor()[1]
#             print(result_comp[result].num_atoms)
#             if ele_res_ratio/ele_ref_ratio > ratio:
#                 ratio = ele_res_ratio/ele_ref_ratio
#         print(ratio,ref_comp[ref_name[reference]].reduced_formula,result_comp[result].reduced_formula)
