#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 12:13:35 2023

@author: bt308570
"""

import sys,os,glob,subprocess,importlib
from pymatgen.core import Structure,Composition

def find_topdir():
    origin = os.getcwd()
    while True:
        if not os.path.isfile(".anchor"):
            os.chdir("..")
        else:
            print("Head directory found.")
            topdir = os.getcwd()
            break
        if os.getcwd() == "/":
            print("No anchor found.")
            sys.exit()
    os.chdir(origin)
    return topdir

sys.stdout.flush()
print("Running script to process completed displacement single-point calculations.")

# locate head directory

workdir = find_topdir()        

# move through the phonopy calculations
phon_dir = workdir+"/VIB"
calcs = [f for f in os.listdir(phon_dir) if os.path.isdir(os.path.join(phon_dir, f))]

for struc in calcs:
    
    # get information on the amount of single point calculations
    os.chdir(phon_dir+"/"+struc)
    structure = Structure.from_file("CONTCAR")
    pattern = 'POSCAR-*'
    files = glob.glob(pattern)
    numbers = [int(os.path.splitext(os.path.basename(file))[0].split('-')[-1]) for file in files]
    n_files = "{:03}".format(max(numbers))

    time = 0.0
    run_phonopy = True
    for i in range(1,max(numbers)+1):
        # check if all calculations are done
        if not os.path.isfile("disp-"+str("{:03}".format(i))+"/done"):
          run_phonopy = False  
        else:
            with open("disp-"+str("{:03}".format(i))+"/OUTCAR", "r") as file:
                for line in file:
                    if "Total C" in line:
                        time = time + float(line.split(":")[1])
          
    # check if the forces were already calculated, run phonopy if not
    if run_phonopy == True and os.path.isfile("FORCE_SETS") == False:
        print("Calculating forces for: "+struc)
        print(subprocess.run(["phonopy -f disp-{001.."+str(n_files)+"}/vasprun.xml"], shell=True, stdout=subprocess.PIPE))
        
    # check if the thermodynamic data was already calculated, run phonopy if not
    #if os.path.isfile("thermal_properties.yaml") == False and
    if os.path.isfile("mesh.conf") == True and os.path.isfile("FORCE_SETS") == True:
        print(subprocess.run(["phonopy mesh.conf -t"], shell=True, stdout=subprocess.PIPE))
        # remove xml files!
    
    # if thermodynamic data is available, extract all necessary information
    if os.path.isfile("thermal_properties.yaml") == True:
        
        temperature = []
        free_energy = []
        t_prop = False
        
        # extracting information
        with open('thermal_properties.yaml','r') as file:
            for line in file:
                if 'thermal_properties:' in line:
                    t_prop = True
                if '- temperature' in line and t_prop == True:
                    temperature.append(line.split()[2])
                if 'free_energy' in line and t_prop == True:
                    free_energy.append(line.split()[1])
        
        # setup everything to process and write results
        composition = Composition(structure.composition)
        Z = composition.get_reduced_composition_and_factor()[1]
        
        spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
        directory_tools = importlib.util.module_from_spec(spec)
        sys.modules["name"] = directory_tools
        spec.loader.exec_module(directory_tools)
        
        res_dir = (workdir+"/RESULTS/"+struc)
        directory_tools.make_directory(res_dir)      
        
        # writing information
        with open(res_dir+'/thermo.txt', 'w') as file:
            file.write("      T        G   G/unit\n")
            file.write("      K       eV  eV/unit\n")
            for i in range(0,len(temperature)):
                t = "{:7.2f}".format(float(temperature[i]))
                g = "{:10.4f}".format(float(free_energy[i])/96.49)
                g_atom = "{:10.4f}".format(float(free_energy[i])/96.49/Z)
                file.write(t+' '+g+" "+g_atom+"\n")
                
        with open(res_dir+'/times.txt', 'w') as file:
            file.write("PHONON   "+"{:10.0f}".format(time)+"\n")
            file.write("per disp "+"{:10.2f}".format(time/max(numbers))+"\n")
            file.write("TOTAL    "+"{:10.0f}".format(time)+"\n")
        
print("\n")