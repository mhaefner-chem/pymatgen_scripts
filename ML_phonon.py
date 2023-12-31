#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 13:31:00 2023

@author: bt308570
"""

from pymatgen.core import Structure,Composition
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Kpoints,Incar,Poscar
from pymatgen.io.vasp.sets import MPScanRelaxSet
from pymatgen.transformations.standard_transformations import SupercellTransformation

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
    print("Running script to create ML-supported phonon calculations.")
    topdir = find_topdir()
    
    T = 400.0
    EDIFF = "-7"
    ACC = True
    ML_WTIFOR = 1.0
    MAXMIX = -1
    alt_settings = ""
    
    # read command line arguments
    n = len(sys.argv)
    for i in range(0,n):
       # reference structure
          if sys.argv[i] == "-t":
              T = int(sys.argv[i+1])
          elif sys.argv[i] == "-e" or sys.argv[i] == "--energy":
              EDIFF = sys.argv[i+1]
          elif sys.argv[i] == "-n" or sys.argv[i] == "--normal":
              ACC = False
          elif sys.argv[i] == "-f" or sys.argv[i] == "--force":
              ML_WTIFOR = float(sys.argv[i+1])
          elif sys.argv[i] == "-m" or sys.argv[i] == "--maxmix":
              MAXMIX = int(sys.argv[i+1])  
              
    alt_settings = alt_settings + "_" + str(T)   
    alt_settings = alt_settings + "_" + EDIFF
    alt_settings = alt_settings + "_" + str(ML_WTIFOR)
    if MAXMIX > 0:
        alt_settings = alt_settings + "_" + str(MAXMIX)
    else:
        alt_settings = alt_settings + "_NOMIX"
    if ACC == True:
        alt_settings = alt_settings + "_ACC"
    else:
        alt_settings = alt_settings + "_NORM"
    
    #     if sys.argv[i] == "-i":
    #         incar = Incar.from_file(sys.argv[i+1])
    #     if sys.argv[i] == "-k":
    #         kpoints = Kpoints.from_file(sys.argv[i+1])
    
    print("Chosen temperature:","{:7.2f}".format(T))
    
    if os.path.isfile("POSCAR"):
        structure = Structure.from_file("POSCAR")
    elif os.path.isfile("CONTCAR"):
        structure = Structure.from_file("CONTCAR")
    else:
        print("Structure file (POSCAR/CONTCAR) missing!")
        sys.exit()
        
    if os.path.isfile("INCAR"):
        incar = Incar.from_file("INCAR")
    else:
        print("INCAR file missing!")
        sys.exit()
        
    if os.path.isfile("KPOINTS"):
        kpoints = Kpoints.from_file("KPOINTS")
    else:
        print("KPOINTS file missing!")
        sys.exit()
    name = os.path.basename(os.getcwd())

# import the directory tools and settings
    spec = importlib.util.spec_from_file_location("directory_tools",topdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)
    
    spec = importlib.util.spec_from_file_location("settings",topdir + "/bin/settings.py")
    settings = importlib.util.module_from_spec(spec)
    sys.modules["name"] = settings
    spec.loader.exec_module(settings)

    tmin = 0.0
    tmax = 1500.0
    tstep = 50.0
    sc_min = 12.0
    learn_min = 4.0
    ML_steps = [250] #[100,250,500] #,1000] #,2500,5000]
    
    time_ml = {}
    for i in ML_steps:
        time_ml[i] = [-1,-1,0]
    
    if os.path.isfile(topdir+"/SETTINGS/settings.txt"):
        with open(topdir+"/SETTINGS/settings.txt",'r') as f:
            lines = f.readlines()
            # name = lines[0]
            for i in range(0,len(lines)):
                if "TMIN" in lines[i]:
                    tmin = float(lines[i].split("=")[1])
                if "TMAX" in lines[i]:
                    tmax = float(lines[i].split("=")[1])
                if "TSTEP" in lines[i]:
                    tstep = float(lines[i].split("=")[1])
                if "SC_MIN" in lines[i]:
                    sc_min = float(lines[i].split("=")[1])


    # take supplied structure

    ML_dir = topdir+"/ML"
    directory_tools.make_directory(ML_dir)

    os.chdir(ML_dir)
    # name = directory_tools.directory_name_constructor(structure)
    workdir = os.path.join(ML_dir,name)
    
    # alternative settings
    # alternate EDIFF
    # alt_settings = "_1E-5"
    # incar["EDIFF"] = 1e-5
    
    # alt_settings = "_FS_1.5"
    # incar["ML_WTIFOR"] = 1.5
    
    
   
    workdir = workdir + alt_settings
    
    directory_tools.make_directory(workdir)
    os.chdir(workdir)
    
    # INCAR manipulation
    
    incar["IBRION"] = 0
    composition = Composition(structure.composition)
    row2 = ["C","N","O","F"]
    row2_elements = []
    for i in row2:
        row2_elements.append(Element(i))
    for i in row2_elements:
        if i in composition.elements:
            incar["POTIM"] = 1.5
            break
        else:
            incar["POTIM"] = 2.0
    
    print(incar["POTIM"])
    
    incar["MDALGO"] = 3
    LG = []
    for i in range(structure.ntypesp):
       LG.append(1.0)
    incar["LANGEVIN_GAMMA"] = LG
    incar["LANGEVIN_GAMMA_L"] = 10
    incar["PMASS"] = 10
    incar["TEBEG"] = T
    incar["ISIF"] = 3
    incar["ML_LMLFF"] = "T"
    incar["ML_WTSIF"] = 2
    incar["RANDOM_SEED"] = [574299384,0,0]
    if MAXMIX > 0:
        incar["MAXMIX"] = MAXMIX
        
    incar["EDIFF"] = "1e"+EDIFF
    incar["ML_WTIFOR"] = ML_WTIFOR

    
    
    unit_struc = structure
    # using a supercell with all lattice parameters > 4 AA each
    lat_param_alt = [0,0,0]
    dim = [1,1,1]
    for i in range(0,3):
        lat_param_alt[i] = structure.lattice.abc[i]
        while lat_param_alt[i] < learn_min:
            lat_param_alt[i] = lat_param_alt[i] + structure.lattice.abc[i]
            dim[i] = dim[i] + 1
    
    st = SupercellTransformation.from_scaling_factors(dim[0],dim[1],dim[2])
    structure = st.apply_transformation(structure)
    
    # write new KPOINTS file with reduced k-points
    reduced_kpts = kpoints.kpts
    for i in range(0,3):
        reduced_kpts[0][i] = math.ceil(kpoints.kpts[0][i]/dim[i])
    kpoints.kpts = reduced_kpts
    kpoints.write_file("KPOINTS")
    
    vasp_input = MPScanRelaxSet(structure)
    vasp_input.write_input(".")
    
    with open("INCAR",'r') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            if "MAGMOM" in lines[i]:
                magmoms = lines[i]
    
    if ACC == True:
        accuracy = 12.0
        incar["NGX"] = math.ceil(structure.lattice.abc[0] * accuracy)
        incar["NGY"] = math.ceil(structure.lattice.abc[1] * accuracy)
        incar["NGZ"] = math.ceil(structure.lattice.abc[2] * accuracy)
    else:
        incar.pop("NGX")
        incar.pop("NGY")
        incar.pop("NGZ")
    
    incar.write_file("INCAR")
    with open("INCAR",'r') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            if "MAGMOM" in lines[i]:
                lines[i] = magmoms
                break
    with open("INCAR",'w') as f:
        for i in range(0,len(lines)):
            f.write(lines[i])
    incar = Incar.from_file("INCAR")
    structure.to("POSCAR",fmt="vasp")
    
    
    poscar = Poscar(structure)
    unit_poscar = Poscar(unit_struc)
    
    for i in ML_steps:
        directory_tools.make_directory(str(i))
        os.chdir(str(i))
        
        directory_tools.make_directory("MD")
        shutil.copyfile(workdir+"/POTCAR","MD/POTCAR") 
        # MD settings
        os.chdir("MD")
        incar["NSW"] = i
        incar["ISYM"] = 0
        incar["ML_ISTART"] = 0
        incar["ML_MODE"] = "train"
        incar.write_file("INCAR_x")
        with open("INCAR_x", "rt") as fin:
            with open("INCAR", "wt") as fout:
                for line in fin:
                    fout.write(line.replace('Train', 'train'))
        kpoints.write_file("KPOINTS")
        poscar.write_file(name+"_"+str(i)+"_MD.inp")
        os.chdir("..")
        
        directory_tools.make_directory("REFIT")
        shutil.copyfile(workdir+"/POTCAR","REFIT/POTCAR")
        os.chdir("REFIT")
        incar["ML_ISTART"] = 1
        incar["ISYM"] = 0
        incar["ML_MODE"] = "refit"
        incar.write_file("INCAR_x")
        with open("INCAR_x", "rt") as fin:
            with open("INCAR", "wt") as fout:
                for line in fin:
                    fout.write(line.replace('Refit', 'refit'))
        kpoints.write_file("KPOINTS")
        poscar.write_file(name+"_"+str(i)+"_REFIT.inp")
        os.chdir("..")
        
        directory_tools.make_directory("VIB")
        shutil.copyfile(workdir+"/POTCAR","VIB/POTCAR")
        os.chdir("VIB")
        incar["ML_ISTART"] = 2
        incar["ML_MODE"] = "run"
        incar["NSW"] = 0
        incar.write_file("INCAR_x")
        with open("INCAR_x", "rt") as fin:
            with open("INCAR", "wt") as fout:
                for line in fin:
                    fout.write(line.replace('Run', 'run'))
        kpoints.write_file("KPOINTS")
        unit_poscar.write_file("POSCAR")
        os.chdir(workdir)
    
    
    # set up ML calculation
    for i in ML_steps:
        os.chdir(str(i))
        if os.path.isfile(os.getcwd()+"/MD/run") or os.path.isfile(os.getcwd()+"/MD/done"):
            print("MD running or done.")
            if os.path.isfile(os.getcwd()+"/MD/done") and os.path.isfile(os.getcwd()+"/MD/OUTCAR"):
                with open(os.getcwd()+"/MD/OUTCAR", "r") as file:
                    for line in file:
                        if "Total C" in line:
                            time_ml[i][0] = float(line.split(":")[1])
        else:
            print(os.getcwd())
            os.chdir("MD")
            incar.write_file("run")
            # run VASP calculation
            print(subprocess.run(["/home/70/bt308570/bin/vsub "+name+"_"+str(i)+"_MD.inp"], shell=True, stdout=subprocess.PIPE))
            os.chdir("..")
        if os.path.isfile(os.getcwd()+"/REFIT/run") or os.path.isfile(os.getcwd()+"/REFIT/done"):
            print("REFIT running or done.")
            if os.path.isfile(os.getcwd()+"/REFIT/done") and os.path.isfile(os.getcwd()+"/REFIT/OUTCAR"):
                with open(os.getcwd()+"/REFIT/OUTCAR", "r") as file:
                    for line in file:
                        if "Total C" in line:
                            time_ml[i][1] = float(line.split(":")[1])
        elif os.path.isfile(os.getcwd()+"/MD/done"):
            if os.path.isfile(os.getcwd()+"/MD/ML_ABN"):
                shutil.copyfile(os.getcwd()+"/MD/ML_ABN",os.getcwd()+"/REFIT/ML_AB")
                os.chdir("REFIT")
                incar.write_file("run")
                # run VASP calculation
                print(subprocess.run(["/home/70/bt308570/bin/vsub "+name+"_"+str(i)+"_REFIT.inp"], shell=True, stdout=subprocess.PIPE))
                os.chdir("..")
        if os.path.isfile(os.getcwd()+"/REFIT/done") and os.path.isfile(os.getcwd()+"/REFIT/ML_FFN"):
            shutil.copyfile(os.getcwd()+"/REFIT/ML_FFN",os.getcwd()+"/VIB/ML_FF")
            os.chdir("VIB")
            
            # create phonopy.conf and mesh.conf 
            
            # using a supercell with all lattice parameters > sc_min AA each
            lat_param_alt = [0,0,0]
            dim = [1,1,1]
            for j in range(0,3):
                lat_param_alt[j] = structure.lattice.abc[j]
                while lat_param_alt[j] < sc_min:
                    lat_param_alt[j] = lat_param_alt[j] + structure.lattice.abc[j]
                    dim[j] = dim[j] + 1
                    
            # write new KPOINTS file with reduced k-points
            reduced_kpts = kpoints.kpts
            for j in range(0,3):
                reduced_kpts[0][j] = math.ceil(kpoints.kpts[0][j]/dim[j])
            kpoints.kpts = reduced_kpts
            kpoints.write_file("KPOINTS")
            
            # write phonopy.conf
            with open('phonopy.conf', 'w') as file:
                file.write("DIM = ")
                for j in range(0,3):
                    file.write(str(dim[j])+" ")
                file.write("\nCREATE_DISPLACEMENTS = .TRUE.\n")
            
            # write mesh.conf
            with open('mesh.conf', 'w') as file:
                file.write("DIM = ")
                for j in range(0,3):
                    file.write(str(dim[j])+" ")
                file.write("\nMP = ")
                # for i in range(0,3):
                    # file.write(str(reduced_kpts[0][i]*4)+" ")
                file.write("15 15 15")
                file.write("\nGAMMA_CENTER = .TRUE.\nTMIN = "+str(tmin)+"\nTMAX = "+str(tmax)+"\nTSTEP = "+str(tstep)+"\n")
            
            # run phonopy to setup the calculations
            if not os.path.isfile(os.getcwd()+"/SPOSCAR"):
                print(subprocess.run(["phonopy phonopy.conf -c POSCAR"], shell=True, stdout=subprocess.PIPE)) 
                
            # perform vasp calculations (adapted from phonopy_vasp_in.sh)
            # identify the required amount of single point calculations
            pattern = 'POSCAR-*'
            files = glob.glob(pattern)
            numbers = []
            numbers = [int(os.path.splitext(os.path.basename(file))[0].split('-')[-1]) for file in files]
            n_files = "{:03}".format(max(numbers))
            
            
            for j in range(1,max(numbers)+1):
                disp_dir = "disp-"+"{:03}".format(j)
                if os.path.isfile(disp_dir+"/run"):
                    print(disp_dir+" was skipped, calculation already exists and is running!")
                elif not os.path.isfile(disp_dir+"/done"):
                    # create necessary folders and copy the relevant files
                    directory_tools.make_directory(disp_dir)
                    shutil.copyfile("KPOINTS", disp_dir+"/KPOINTS")
                    shutil.copyfile("POTCAR", disp_dir+"/POTCAR")
                    shutil.copyfile("ML_FF", disp_dir+"/ML_FF")
                    poscar_name = name+"-"+"{:03}".format(j)+".inp"
                    shutil.copyfile("POSCAR-"+"{:03}".format(j), disp_dir+"/"+poscar_name)
                    shutil.copyfile("KPOINTS", disp_dir+"/KPOINTS")
                    # dirty fix for magmom
                    with open("INCAR", "rt") as fin:
                        with open("INCAR_super", "wt") as fout:
                            for line in fin:
                                if "MAGMOM" in line:
                                    line = "MAGMOM = 10000*0.6\n"
                                fout.write(line)
                    shutil.copyfile("INCAR_super", disp_dir+"/INCAR")
                    
                    tmp_struc = Structure.from_file("POSCAR-"+"{:03}".format(j))
                    
                    os.chdir(disp_dir)


                    shutil.copyfile(poscar_name,"POSCAR")
                    # run VASP calculation
                    print(subprocess.run(["/home/70/bt308570/vasp.6.4.1_neb/bin/vasp_std"], shell=True, stdout=subprocess.PIPE))
                    incar.write_file("done")
                    os.chdir("..")
                elif os.path.isfile(disp_dir+"/done"):
                    with open(disp_dir+"/OUTCAR", "r") as file:
                        for line in file:
                            if "Total C" in line:
                                time_ml[i][2] = time_ml[i][2] + float(line.split(":")[1])
                    
            run_phonopy = True
            for j in range(1,max(numbers)+1):
                # check if all calculations are done
                if not os.path.isfile("disp-"+str("{:03}".format(j))+"/done"):
                  run_phonopy = False  
                  
            # check if the forces were already calculated, run phonopy if not
            if run_phonopy == True and os.path.isfile("FORCE_SETS") == False:
                print("Calculating forces for: "+name)
                print(subprocess.run(["phonopy -f disp-{001.."+str(n_files)+"}/vasprun.xml"], shell=True, stdout=subprocess.PIPE))
                
            # check if the thermodynamic data was already calculated, run phonopy if not
            #if os.path.isfile("thermal_properties.yaml") == False and
            if os.path.isfile("mesh.conf") == True and os.path.isfile("FORCE_SETS") == True and os.path.isfile("thermal_properties.yaml") == False:
                print(subprocess.run(["phonopy mesh.conf -t"], shell=True, stdout=subprocess.PIPE))
                
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
                
                res_dir = (topdir+"/RESULTS/"+name+"_ML_"+str(i)+alt_settings)
                directory_tools.make_directory(res_dir)      
                
                shutil.copy(topdir+"/RESULTS/"+name+"/energy.txt",res_dir+"/energy.txt")
                
                # writing information
                with open(res_dir+'/thermo.txt', 'w') as file:
                    file.write("      T        G   G/unit\n")
                    file.write("      K       eV  eV/unit\n")
                    for temps in range(0,len(temperature)):
                        t = "{:7.2f}".format(float(temperature[temps]))
                        g = "{:10.4f}".format(float(free_energy[temps])/96.49)
                        g_atom = "{:10.4f}".format(float(free_energy[temps])/96.49/Z)
                        file.write(t+' '+g+" "+g_atom+"\n")
                        
        
                        
                with open(res_dir+'/times.txt', 'w') as file:
                    file.write("MD       "+"{:10.0f}".format(time_ml[i][0])+"\n")
                    file.write("REFIT    "+"{:10.0f}".format(time_ml[i][1])+"\n")
                    file.write("FF       "+"{:10.0f}".format(time_ml[i][2])+"\n")
                    file.write("per disp "+"{:10.2f}".format(time_ml[i][2]/max(numbers))+"\n")
                    total = time_ml[i][0] + time_ml[i][1] + time_ml[i][2]
                    file.write("TOTAL    "+"{:10.0f}".format(total)+"\n")
        
        os.chdir(workdir)
    
    print("Ta-dah!")
    

    
    # make subdirectory named after structure - maybe make supercell if lattice parameter < 6 AA
    # make ML directories (MD, Refit, VIB)
    # do vibrational computations once refit was carried out

    sys.exit()
