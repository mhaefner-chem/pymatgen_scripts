#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:34:18 2023

@author: bt308570
"""

import sys,os,shutil,inspect
from pathlib import Path

# this function creates a directory
def make_directory(dir):
    import os
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        print('The directory '+os.path.basename(dir)+' already exists and was skipped.')
    return

# this function returns some information on the program and how to run it
def program_help():
    print("**********\nThis script sets up the folder environment for a new phase analysis. It requires a name and all elements involved in the calculations.\n")
    print("A custom name can be supplied with '-n NAME' or '--name NAME'.")
    print("The elements are supplied with '-e #N LIST_OF_ELEMENTS' as in '-e 3 Na Cl Mg'.\n")
    print("This help dialogue is displayed with either '-h' or '--help'.")
    print("**********")
    sys.exit()

if __name__ == '__main__':
    sys.stdout.flush()
    element = []
    
    # read in all arguments
    n = len(sys.argv)
    for i in range(0,n):
        if sys.argv[i] == "-n" or sys.argv[i] == "--name":
            name = sys.argv[i+1]
        elif sys.argv[i] == "-e" or sys.argv[i] == "--element":
            n_element = int(sys.argv[i+1])
            for j in range(0,n_element):
                element.append(sys.argv[i+j+2])
        elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
            program_help()

    # check if necessary information was supplied, stop script otherwise
    if not "name" in globals():
        print("No name supplied with '-n'. Stopping script.")
        sys.exit()
    if not "n_element" in globals():
        print("No elements supplied with '-e #elements H He Li...'. Stopping script.")
        sys.exit()            
    
    # setup folder network   
    
    make_directory(name)
    os.chdir(name)
    workdir = os.getcwd()
    
    make_directory(workdir+"/CIF/EXP")
    make_directory(workdir+"/CIF/REF")
    # make_directory(workdir+"/PREP/CONV")
    # make_directory(workdir+"/PREP/SUPER")
    # make_directory(workdir+"/PREP/PREOPT")
    make_directory(workdir+"/OPT")
    make_directory(workdir+"/VIB")
    make_directory(workdir+"/RESULTS")
    make_directory(workdir+"/SETTINGS")
    make_directory(workdir+"/bin")
    
    # copies relevant scripts from its origin
    origin = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    files = ["directory_tools.py","settings.py","converge_cifs.py","phonocalc.py","phonoresult.py","result_collector.py","env_run.py","optimizer.py","supercell_preopt.py"]
    for file in files:
        shutil.copyfile(origin+"/"+file,workdir+"/bin/"+file)
    
    
    # specifies an anchor for other programs
    Path(workdir+"/.anchor").touch()
    #Path(workdir+"/SETTINGS/settings.txt").touch()
    
    # writes the specified elements to a file in SETTINGS, used later on for the SHIFT module
    with open('SETTINGS/element.txt', 'w') as file:
        for i in element:
         file.write(i+"\n")
    
sys.exit()
