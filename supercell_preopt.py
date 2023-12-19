#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:26:54 2023

@author: bt308570
"""

from pymatgen.core import Structure

import sys,os,shutil
from glob import glob
import importlib.util
import subprocess,signal,warnings
from pathlib import Path

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
    
    file_dir = os.getcwd()
    
    selec = "100"
    super = 0
    mode = "vasp"
    n = len(sys.argv)
    for i in range(0,n):
        # specify the cells and supercell size from cif
        if sys.argv[i] == "-c":
            print("Reading in the file "+sys.argv[i+1]+".")
            structure = Structure.from_file(sys.argv[i+1])
            cif_file = sys.argv[i+1]
            super = 1
        elif sys.argv[i] == "-r" or sys.argv[i] == "--reference":
            print("Reading in the file "+sys.argv[i+1]+".")
            reference = Structure.from_file(sys.argv[i+1])
            ref_file = sys.argv[i+1]
        elif sys.argv[i] == "-d" or sys.argv[i] == "--dopant":
            print("Reading in the file "+sys.argv[i+1]+".")
            dopant = Structure.from_file(sys.argv[i+1])
            dope_file = sys.argv[i+1]
        elif sys.argv[i] == "-s" or sys.argv[i] == "--supercell":
            sc_size = [int(x) for x in sys.argv[i+1].split("x")]
            sc_str = sys.argv[i+1]
        elif sys.argv[i] == "-n" or sys.argv[i] == "--name":
            id_tag = sys.argv[i+1]
        elif sys.argv[i] == "-m" or sys.argv[i] == "--mode":
            mode = sys.argv[i+1]
        elif sys.argv[i] == "-q" or sys.argv[i] == "--coulomb":
            selec = sys.argv[i+1]
        elif sys.argv[i] == "-p" or sys.argv[i] == "--potential":
            mlp_file = sys.argv[i+1]
        
    
    warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
    warnings.filterwarnings("ignore", message="POTCAR data with symbol") # suppress the POTCAR warning, works regardless
    
    
    
    
    print("Running script to carry out supercell calculations.")
    # locate head directory
    workdir = find_topdir()
    os.chdir(workdir)
    
    # import the directory tools 
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)   
    
    directory_tools.make_directory("SUPERCELL")
    os.chdir("SUPERCELL")
    directory_tools.make_directory("PREOPT")
    directory_tools.make_directory("OPT")
    os.chdir("PREOPT")
    
    
    # set up the calculation for the supercell-generated structures
    if super == 1:
        name_prefix = directory_tools.directory_name_constructor(structure)
        name_suffix = directory_tools.directory_name_constructor(reference)
        
        name = name_prefix+"_"+sc_str+"_"+id_tag+"_base_"+name_suffix
        directory_tools.make_directory(name)
        os.chdir(name)
        
        workdir = os.getcwd()
        
        shutil.copyfile(os.path.join(file_dir,cif_file),cif_file)
        shutil.copyfile(os.path.join(file_dir,mlp_file),mlp_file)
        
        # run supercell to generate all relevant structures
        if not os.path.isfile("done"):
            process = subprocess.Popen(['supercell',"-i",cif_file,"-s",sc_str,"-n","l"+selec,"-q","-m"],
                         stdin=subprocess.DEVNULL,
                         start_new_session=False)
            process.wait()
            shutil.copyfile(cif_file,"done")
    
        
        files = glob("supercell*.cif")
        
        # shutil.copyfile(os.path.join(file_dir,ref_file),"ref.cif")
        # shutil.copyfile(os.path.join(file_dir,dope_file),"dope.cif")
        
        counter = 0
        for i in files:
            base_name = i.split(".")[0]
            if mode == "szp":
                base_name = base_name + "_gpaw_szp"
            elif mode == "dzp":
                base_name = base_name + "_gpaw_dzp"
            elif mode == "fast":
                  base_name = base_name + "_fast"
            elif mode == "slow":
                  base_name = base_name + "_slow"
            elif mode == "ml":
                  base_name = base_name + "_ml"
                   
            directory_tools.make_directory(base_name)
            shutil.copyfile(i,os.path.join(base_name,i))
            if mode == "ml":
                shutil.copyfile(mlp_file,os.path.join(base_name,mlp_file))
            os.chdir(base_name)
            
            if not os.path.isfile("INCAR"):
                print("Setting up input!")
                counter += 1
                print(counter)
                if mode == "fast":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--fast","--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                elif mode == "slow":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                elif mode == "ml":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--ionly","--ml"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                    shutil.copyfile("INCAR","INCAR_x")
                    with open("INCAR_x", "rt") as fin:
                        with open("INCAR", "wt") as fout:
                            for line in fin:
                                fout.write(line.replace('Run', 'run'))
                else:
                    
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"-p","--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
            else:
                print("Already set up.")
            
            os.chdir(workdir)
        
    # set up the calculations for the reference structure
    elif super == 0:
        os.chdir(os.path.join(find_topdir(),"SUPERCELL","PREOPT"))
        name = directory_tools.directory_name_constructor(reference)
        directory_tools.make_directory(name)
        os.chdir(name)
        shutil.copyfile(os.path.join(file_dir,ref_file),ref_file)
        process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",ref_file,"-n",name,"-p"],
                     stdin=subprocess.DEVNULL,
                     stdout=open(name+'.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=True,
                     preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    
    # set up the calculations for the dopant structure
    
    # os.chdir(os.path.join(find_topdir(),"SUPERCELL","PREOPT"))
    # name = directory_tools.directory_name_constructor(dopant)
    # directory_tools.make_directory(name)
    # os.chdir(name)
    # shutil.copyfile(os.path.join(file_dir,dope_file),dope_file)
    # process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",dope_file,"-n",name,"-p"],
    #              stdin=subprocess.DEVNULL,
    #              stdout=open(name+'.out', 'w'),
    #              stderr=subprocess.STDOUT,
    #              start_new_session=True,
    #              preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    
    
    sys.exit()
    
