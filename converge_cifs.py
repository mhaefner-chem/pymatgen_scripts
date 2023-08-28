#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:26:54 2023

@author: bt308570
"""

from pymatgen.core import Structure,Composition
from pymatgen.core.periodic_table import Species
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun

import sys,os,shutil
from glob import glob
import importlib.util
import subprocess,signal,warnings

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

# beginning of main program
if __name__ == '__main__':
    sys.stdout.flush()
    warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
    warnings.filterwarnings("ignore", message="POTCAR data with symbol") # suppress the POTCAR warning, works regardless

    print("Running script to converge existing cif-files.")
    # locate head directory
    workdir = find_topdir()
    os.chdir(workdir)
    
    # find cif-files
    exp_dir = workdir+"/CIF/EXP"
    ref_dir = workdir+"/CIF/REF"
    element_dir = workdir+"/SETTINGS"
    exp_files = glob(exp_dir + "/*cif")
    sc_files = glob(workdir + "/SUPERCELL/BEST/*vasp")
    ref_files = glob(ref_dir + "/*.cif")
    structures = {}
    elements =  {}
    oxi_states = {}
    name = {}
    name_suffixes = {}
    
    # import the directory tools 
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)
    
    # read in the elements for the swap
    with open(element_dir+'/element.txt','r') as file:
        for line in file:
            species = Species.from_string(line)
            ele_name = species.element.symbol
            elements[ele_name] = 0
            oxi_states[ele_name] = species.oxi_state
    
    # read in cif files for the reference structures
    for ref_file in ref_files:
        print(ref_file)
        structure = Structure.from_file(ref_file)
        composition = Composition(structure.composition)
        name_suffix = composition.reduced_formula
        name_suffix = name_suffix.replace('(','_')
        name_suffix = name_suffix.replace(')','_')
        key = os.path.basename(os.path.splitext(ref_file)[0])
        structures[key] = structure
        name_suffixes[key] = name_suffix
        print("Working on: "+key)

    
    # substitute the elements
    for structure in structures:
        print(structure,", substituting")
        for element in elements:
            elements[element] = 0        
        # count the occurrences of the different elements
        for site in structures[structure].sites:
            if site.specie.symbol in elements:
                elements[site.specie.symbol] = elements[site.specie.symbol] + 1
        
        # identify the missing elements in the structure
        switch_in = []
        for i in elements:
            if elements[i] == 0:
                switch_in.append(i)

              
        # swap the divergent element for the missing element
        i = -1
        for site in structures[structure].sites:
            i = i + 1
            if not site.specie.symbol in elements:
                for switcher in switch_in:
                    if abs(site.specie.oxi_state-oxi_states[switcher]) < 1e-3:
                        structures[structure].replace(i,switcher)
                    else:
                        print("Oxi. state mismatch.")
                
        name_prefix = directory_tools.directory_name_constructor(structures[structure])
        name[structure] = name_prefix+"_base_"+name_suffixes[structure]
        
        
    # retrieve structures from the exp-files
    for exp_file in exp_files:
        structure = Structure.from_file(exp_file)
        key = os.path.basename(os.path.splitext(exp_file)[0])
        structures[key] = structure
        print("Working on: "+key)
        name[key] = directory_tools.directory_name_constructor(structure)
    
    for sc_file in sc_files:
        structure = Structure.from_file(sc_file)
        key = os.path.basename(os.path.splitext(sc_file)[0])
        structures[key] = structure
        print("Working on: "+key)
        name[key] = os.path.basename(os.path.splitext(sc_file)[0])
for s in structures:
    print(s)

for structure in structures:
    # assign naming
    struc_dir = (workdir+"/CONV/"+name[structure])
    opt_dir = (workdir+"/OPT/"+name[structure])
    res_dir = (workdir+"/RESULTS/"+name[structure])
    directory_tools.make_directory(struc_dir)
    os.chdir(struc_dir)
    
    # if there is no optimization running, start a new one
    if os.path.isfile('CONVOPT/run') or os.path.isfile('run'):
        print("Optimization " + name[structure] + " is already in the works!")
            
            
    elif os.path.isfile('CONVOPT/done'):
        print("Optimization " + name[structure] + " is done!")
        directory_tools.make_directory(opt_dir)
        files = ["INCAR","KPOINTS","CONTCAR"]
        for file in files:
            if os.path.isfile(struc_dir+"/CONVOPT/"+file):
                shutil.copyfile(struc_dir+"/CONVOPT/"+file,opt_dir+"/"+file)
            else:
                print("The file "+file+" does not exist in directory "+struc_dir+"/CONVOPT!")
        
        # collect all results in RESULTS - energy, energy/unit cell, energy/atom
        xml = Vasprun("vasprun.xml")
        e = "{:10.4f}".format(xml.final_energy)
        e_atom = "{:10.4f}".format(xml.final_energy/structures[structure].num_sites)
                
        composition = Composition(structures[structure].composition)
        Z = composition.get_reduced_composition_and_factor()[1]
        e_Z = "{:10.4f}".format(xml.final_energy/Z)
        
        if not os.path.isdir(res_dir):
            directory_tools.make_directory(res_dir)
        with open(res_dir+'/energy.txt', 'w') as file:
            file.write("                    eV\n")
            file.write("Energy      "+e+"\n")
            file.write("Energy/unit "+e_Z+"\n")
            file.write("Energy/atom "+e_atom+"\n")
        # with open(res_dir+'/charges.txt', 'w') as file:
        #     file.write(site.specie.symbol,"{:4.1f}".format(site.specie.oxi_state)+"\n")
            

        
    else:
        print("Setting up optimization of "+name[structure]+".")
        poscar = Poscar(structures[structure])
        poscar.write_file("POSCAR")
        
        # run optimizer.py to run the VASP calculation
        subprocess.Popen(['python3','-u', workdir+'/bin/optimizer.py','--fullopt','-c','POSCAR','-n',name[structure]],
                     stdin=subprocess.DEVNULL,
                     stdout=open(name[structure]+'.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=True,
                     preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    print("\n")
    os.chdir(workdir)
