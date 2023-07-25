#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:26:54 2023

@author: bt308570
"""
from pymatgen.core import Structure,Composition
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun

import sys,os,shutil
from glob import glob
import importlib.util
import subprocess,signal

# beginning of main program
if __name__ == '__main__':
    sys.stdout.flush()

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
    # find CONTCAR files
    opt_dir = workdir+"/OPT"
    opt_files = glob(opt_dir + "/*/vasprun.xml")
    structures = {}

for opt_file in opt_files:
    structure = Structure.from_file(opt_file)
    key = os.path.basename(os.path.splitext(opt_file)[0])
    structures[key] = structure
    print(key)
    
    xml = Vasprun(opt_file)
        
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)


    name = directory_tools.directory_name_constructor(structure)
    print(name)
    res_dir = (workdir+"/RESULTS/"+name)
    
    # collect results in RESULTS - e.g., the energy and energy/atom
    e = "{:10.4f}".format(xml.final_energy)
    e_atom = "{:10.4f}".format(xml.final_energy/structure.num_sites)
    
    composition = Composition(structure.composition)
    Z = composition.get_reduced_composition_and_factor()[1]
    e_Z = "{:10.4f}".format(xml.final_energy/Z)
    
    if not os.path.isdir(res_dir):
        directory_tools.make_directory(res_dir)
    with open(res_dir+'/energy.txt', 'w') as file:
        file.write("                    eV\n")
        file.write("Energy      "+e+"\n")
        file.write("Energy/unit "+e_Z+"\n")
        file.write("Energy/atom "+e_atom+"\n")
        
    os.chdir(workdir)


