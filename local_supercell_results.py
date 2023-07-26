#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:21:39 2023

@author: bt308570
"""

# this program collects all results and compiles them into one file for each structure

import sys,os,importlib,glob,warnings
from pymatgen.core import Structure,Composition
from pymatgen.io.vasp.outputs import Vasprun
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import statistics


def energy(outcar):
    last_match = None

    
    with open(outcar, 'r') as file:
        for line in file:
            if "free  e" in line:
                last_match = line
            if "Total CPU" in line:
                timing = float(line.split(":")[1])
        
    last_match = last_match.split("=")[1]
    last_match = last_match.split("eV")[0]
    return float(last_match),timing


# beginning of main program

#######################################

e_cut = 0.6

#######################################

sys.stdout.flush()
warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
warnings.filterwarnings("ignore", message="your POTCAR is corrupted")
print("Collecting Results.")

# read command line arguments
n = len(sys.argv)
e_ref = 0
for i in range(0,n):
    # reference energy
    if sys.argv[i] == "-r":
        print("Reference energy: "+sys.argv[i+1]+".")
        e_ref = float(sys.argv[i+1])
    if sys.argv[i] == "-c":
        e_cut = float(sys.argv[i+1])

preopt_dir = os.getcwd()
while True:
    if not os.path.isfile(".anchor"):
        os.chdir("..")
    else:
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
os.chdir(preopt_dir)



###############################################################
# try to remove

res_energy = {}
res_struc = {}
res_comp = {}
results = []

config_energy = {}
timings = {}
config_sp_energy = {}
ml_energy = {}
config_part_energy = {}
ref_energy = {}

###############################################################

# get all configurations
configs = [config for config in os.listdir(os.path.join(preopt_dir)) if os.path.isdir(os.path.join(preopt_dir,config))]
structures = {}

#specify the different calculation types
e_types = ["fullAI","Q","SP","ML","partOpt"] # coulomb, sp, ML, partopt
e_min = [1] * len(e_types)

e_all = {}
for config in configs:
    e_all[config] = [1] * len(e_types)


# read in coulomb reference energies
with open("supercell_coulomb_energy_l.txt","r") as file:
    for line in file:
        temp = line.split("\t")[1]
        index = e_types.index("Q")
        e_all[config][index] = float(temp.split(" ")[0])
        e_min[index] = min(e_all[config][index],e_min[index])


# get #atoms, composition, and Z from one of the structures (identical among all)
structure = Structure.from_file(config+"/POSCAR")
composition = Composition(structure.composition)
Z = composition.get_reduced_composition_and_factor()[1]
atoms = structure.num_sites

for config in configs:
    
    os.chdir(config)
    if os.path.isfile("done") and os.path.isfile("OUTCAR"):
        
        # get full ab-initio results
        
        index = e_types.index("fullAI")
        final_energy,timing = energy("OUTCAR")
        timings[config] = timing
        
        e_all[config][index] = final_energy/Z
        e_min[index] = min(e_all[config][index],e_min[index])
        
        # energy after single-point calculation
        
        index = e_types.index("SP")
        e_all[config][index] = energy("SCFCONV/OUTCAR")[0]/Z
        e_min[index] = min(e_all[config][index],e_min[index])
        
        
        # energy after at max 10 optimization steps
            
        index = e_types.index("partOpt")
        
        i = 0
        with open("OUTCAR", 'r') as file:
            for line in file:
                if "free  e" in line and i < 10:
                    i+=1
                    temp = line
            
        temp = temp.split("=")[1]
        temp = temp.split("eV")[0]
        e_all[config][index] = float(temp)/Z
        e_min[index] = min(e_all[config][index],e_min[index])
        
        # read in ML if it exists
        if os.path.isdir("ML"):
            index = e_types.index("ML")
            e_all[config][index] = energy("ML/OUTCAR")[0]/Z
        else:
            e_all[config][index] = 100
        e_min[index] = min(e_all[config][index],e_min[index])
    
    # elif for GPAW results
    # read out results and timings in gpaw_e and gpaw_t
        
    else:
        e_all[config] = 1000.0
        if os.path.isfile("SCFCONV/done"):
            if os.path.isfile("SCFCONV/vasprun.xml"):
                xml = Vasprun("SCFCONV/vasprun.xml")
            else:
                print("Missing xml for",config)
            # if xml.converged_electronic == False:
            #     print("The SCF convergence failed!")                
        else:
            print(config,"Config not done!")
    
    
    os.chdir(os.path.join(preopt_dir))

# DO REGRESSION STUFF

# regression for coulomb prediction vs. optimized
x_cou = []
x_uopt = []
x_ml = []
x_part = []
y = []
r_ml = 100


for key, values in e_all.items():

    y.append(float(values[0])/y_min)
    x_cou.append(float(values[1])/x_cou_min)
    x_uopt.append(float(values[2])/x_uopt_min)
    x_ml.append(float(values[3])/x_ml_min)
    x_part.append(float(values[4])/x_part_min)
    
x_lin = np.linspace(0, 2, 100)

fig, ax = plt.subplots()

size = 25
linewidth = 0.5
    
if x_cou[0] != x_cou[1]:
    slope_cou, intercept_cou, r, p, std_err = stats.linregress(x_cou,y)
if x_uopt[0] != x_uopt[1]:
    slope_uopt, intercept_uopt, r_opt, p, std_err = stats.linregress(x_uopt,y)
if x_ml[0] != x_ml[1]:
    slope_ml, intercept_ml, r_ml, p, std_err = stats.linregress(x_ml,y)
    y_ml_lin = slope_ml * x_lin + intercept_ml
    plt.scatter(x_ml, y,label="ML-FF,r² = "+"{:6.3f}".format(r_ml**2),edgecolors="#000000",s=size,linewidth=linewidth)
    plt.plot(x_lin, y_ml_lin)
if x_part[0] != x_part[1]:
    slope_part, intercept_part, r_part, p, std_err = stats.linregress(x_part,y)



# def myfunc(x):
#     return slope * x + intercept

# mymodel = list(map(myfunc, x))



y_cou_lin = slope_cou * x_lin + intercept_cou
y_uopt_lin = slope_uopt * x_lin + intercept_uopt
y_part_lin = slope_part * x_lin + intercept_part

plt.scatter(x_cou, y,label="Coulombic only,r² = "+"{:6.3f}".format(r**2),edgecolors="#000000",s=size,linewidth=linewidth)
plt.scatter(x_uopt, y,label="Unoptimized energy,r² = "+"{:6.3f}".format(r_opt**2),edgecolors="#000000",s=size,linewidth=linewidth)
plt.scatter(x_part, y,label="10 Optsteps,r² = "+"{:6.3f}".format(r_part**2),edgecolors="#000000",s=size,linewidth=linewidth)
plt.plot(x_lin, y_cou_lin)
plt.plot(x_lin, y_uopt_lin)

plt.plot(x_lin, y_part_lin)
plt.plot(x_lin, x_lin)

plt.ylim(ymin=0.95, ymax=1.01)
plt.xlim(xmin=0.95, xmax=1.01)
plt.legend()
ax.set_xlabel('rel. Energy/Predictors')
ax.set_ylabel('rel. Energy/PBEsol')
 
plt.show() 


# results
print("r² cou=","{:6.3f}".format(r**2),"Data points:",len(configs))
print("r² uopt=","{:6.3f}".format(r_opt**2),"Data points:",len(configs))
if r_ml < 50:
    print("r² ML=","{:6.3f}".format(r_ml**2),"Data points:",len(configs))

min_energy = 0.0
for key,value in config_energy.items():
    if value < min_energy:
        min_energy = value
        
form = "{:7.3f}"


if e_ref < -1:  
    print("{:24}".format("Configuration"),"{:^15}".format("E-E_min"),"{:^15}".format("E-E_ref"),"{:^7}".format("E"),"{:^7}".format("Del_opt"))
    print("{:24}".format(""),"   cell","   atom", "   cell", "   atom", "   cell")
    for key,value in config_energy.items():
        if value-e_ref > e_cut:
            continue
        delta_e_min_u = form.format(value-min_energy)
        delta_e_min_a = form.format((value-min_energy)*Z/atoms)
        delta_e_ref_u = form.format(value-e_ref)
        delta_e_ref_a = form.format((value-e_ref)*Z/atoms)
        delta_e_opt_u = value-config_sp_energy[key]
        print("{:24}".format(key),delta_e_min_u,delta_e_min_a,delta_e_ref_u,delta_e_ref_a,form.format(value),form.format(delta_e_opt_u),"{:10.2f}".format(timings[key]))
        res_struc[key].to(key+".vasp",fmt="POSCAR")
else:
    print("{:24}".format("Configuration"),"{:^15}".format("E-E_min"),"{:^7}".format("E_tot"),"{:^10}".format("CPU time"))
    print("{:24}".format(""),"   cell","   atom")
    for key,value in config_energy.items():
        delta_e_min_u = form.format(value-min_energy)
        delta_e_min_a = form.format((value-min_energy)*Z/atoms)
        delta_e_ref_u = form.format(value-e_ref)
        delta_e_ref_a = form.format((value-e_ref)*Z/atoms)
        print("{:24}".format(key),delta_e_min_u,delta_e_min_a,"{:6.3f}".format(value),"{:10.2f}".format(timings[key]))
            

# print(min(y))
# print(statistics.mean(y)-min(y))
