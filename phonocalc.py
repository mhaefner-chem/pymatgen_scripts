#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:26:54 2023

@author: bt308570
"""

from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Kpoints,Incar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPScanRelaxSet

import sys,os,shutil,math,glob
import importlib.util
import subprocess

# locate head directory from a subdirectory
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
    print("Running script to create phonon calculations.")
    workdir = find_topdir()
    os.chdir(workdir)

# import the directory tools and settings
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)
    
    spec = importlib.util.spec_from_file_location("settings",workdir + "/bin/settings.py")
    settings = importlib.util.module_from_spec(spec)
    sys.modules["name"] = settings
    spec.loader.exec_module(settings)

    tmin = 0.0
    tmax = 1500.0
    tstep = 50.0
    sc_min = 12.0
    
    if os.path.isfile(workdir+"/SETTINGS/settings.txt"):
        print(1)
        with open(workdir+"/SETTINGS/settings.txt",'r') as f:
            lines = f.readlines()
            name = lines[0]
            for i in range(0,len(lines)):
                if "TMIN" in lines[i]:
                    tmin = float(lines[i].split("=")[1])
                if "TMAX" in lines[i]:
                    tmax = float(lines[i].split("=")[1])
                if "TSTEP" in lines[i]:
                    tstep = float(lines[i].split("=")[1])
                if "SC_MIN" in lines[i]:
                    sc_min = float(lines[i].split("=")[1])
                #     for j in range(i+1,len(lines)):
                #         thermo_g.append(lines[j].split())
                # elif "Delta_E_0" in lines[i]:
                #     pd_struc_deltaE[name] = lines[i].split()[1]
                # elif "SG" in lines[i]:
                #     # determine the reference states for each structure
                #     ref_states[lines[i].split()[0]] = float(lines[i].split()[1])
                #     if not lines[i].split()[0] in all_ref_states:
                #         all_ref_states.append(lines[i].split()[0])

    # find converged optimizations

    opt_dir = workdir+"/OPT"
    phon_dir = workdir+"/VIB"
    opts = [f for f in os.listdir(opt_dir) if os.path.isdir(os.path.join(opt_dir, f))]

    for struc in opts:
        print(struc)
        
        # copy the directory with the optimized structure
        if not os.path.isdir(phon_dir+"/"+struc):
            shutil.copytree(opt_dir+"/"+struc,phon_dir+"/"+struc)     
        os.chdir(phon_dir+"/"+struc)
        
        # read the structure and retrieve the KPOINTS file from the optimization
        if os.path.isfile("CONTCAR"):
            structure = Structure.from_file("CONTCAR")
        if not os.path.isfile("KPOINTS") and os.path.isfile(opt_dir+"/"+struc+"/KPOINTS"):
            shutil.copyfile(opt_dir+"/"+struc+"/KPOINTS","KPOINTS") 
        elif os.path.isfile("KPOINTS") and os.path.isfile(opt_dir+"/"+struc+"/KPOINTS"):
            os.remove("KPOINTS")
            shutil.copyfile(opt_dir+"/"+struc+"/KPOINTS","KPOINTS") 

        
        # create phonopy.conf and mesh.conf 
        if os.path.isfile("KPOINTS"):
            kpoints = Kpoints.from_file("KPOINTS")
        elif os.path.isfile("vasprun.xml"):
            xml = Vasprun("vasprun.xml")
            kpoints = xml.kpoints
        
        # using a supercell with all lattice parameters > 10 AA each
        lat_param_alt = [0,0,0]
        dim = [1,1,1]
        for i in range(0,3):
            lat_param_alt[i] = structure.lattice.abc[i]
            while lat_param_alt[i] < sc_min:
                lat_param_alt[i] = lat_param_alt[i] + structure.lattice.abc[i]
                dim[i] = dim[i] + 1
                
        # write new KPOINTS file with reduced k-points
        reduced_kpts = kpoints.kpts
        for i in range(0,3):
            reduced_kpts[0][i] = math.ceil(kpoints.kpts[0][i]/dim[i])
        kpoints.kpts = reduced_kpts
        kpoints.write_file("KPOINTS")
        
        # write phonopy.conf
        with open('phonopy.conf', 'w') as file:
            file.write("DIM = ")
            for i in range(0,3):
                file.write(str(dim[i])+" ")
            file.write("\nCREATE_DISPLACEMENTS = .TRUE.\n")
        
        # write mesh.conf
        with open('mesh.conf', 'w') as file:
            file.write("DIM = ")
            for i in range(0,3):
                file.write(str(dim[i])+" ")
            file.write("\nMP = ")
            # for i in range(0,3):
                # file.write(str(reduced_kpts[0][i]*4)+" ")
            file.write("15 15 15")
            file.write("\nGAMMA_CENTER = .TRUE.\nTMIN = "+str(tmin)+"\nTMAX = "+str(tmax)+"\nTSTEP = "+str(tstep)+"\n")
        
        # run phonopy to setup the calculations
        if not os.path.isfile(phon_dir+"/"+struc+"/SPOSCAR"):
            print(subprocess.run(["phonopy phonopy.conf -c CONTCAR"], shell=True, stdout=subprocess.PIPE)) 
            os.chdir(phon_dir)
            
        # perform vasp calculations (adapted from phonopy_vasp_in.sh)
        # identify the required amount of single point calculations
        os.chdir(phon_dir+"/"+struc)
        pattern = 'POSCAR-*'
        files = glob.glob(pattern)
        numbers = []
        numbers = [int(os.path.splitext(os.path.basename(file))[0].split('-')[-1]) for file in files]
        n_files = "{:03}".format(max(numbers))
        
        
        for i in range(1,max(numbers)+1):
            disp_dir = "disp-"+"{:03}".format(i)
            if os.path.isfile(disp_dir+"/run"):
                print(struc+"/"+disp_dir+" was skipped, calculation already exists and is running!")
            elif not os.path.isfile(disp_dir+"/done"):
                # create necessary folders and copy the relevant files
                directory_tools.make_directory(disp_dir)
                shutil.copyfile("KPOINTS", disp_dir+"/KPOINTS")
                poscar_name = struc+"-"+"{:03}".format(i)+".inp"
                shutil.copyfile("POSCAR-"+"{:03}".format(i), disp_dir+"/"+poscar_name)
                
                tmp_struc = Structure.from_file("POSCAR-"+"{:03}".format(i))
                
                os.chdir(disp_dir)
                
                # adapt INCAR file, generate POTCAR
                vasp_input = MPScanRelaxSet(tmp_struc)
                vasp_input.user_incar_settings = settings.incar_settings("phonopy_sps")
                accuracy = 12.0
                vasp_input.user_incar_settings["NGX"] = math.ceil(tmp_struc.lattice.abc[0] * accuracy)
                vasp_input.user_incar_settings["NGY"] = math.ceil(tmp_struc.lattice.abc[1] * accuracy)
                vasp_input.user_incar_settings["NGZ"] = math.ceil(tmp_struc.lattice.abc[2] * accuracy)
                vasp_input.write_input(".")
                incar = Incar.from_file("INCAR")
                if "GGA" in incar and "METAGGA" in incar:
                    incar.pop("METAGGA")
                incar.write_file("INCAR")
                incar.write_file("run")
                
                # run VASP calculation
                print(subprocess.run(["/home/70/bt308570/bin/vsub_py "+poscar_name], shell=True, stdout=subprocess.PIPE))
                os.chdir(phon_dir+"/"+struc)
print("\n")