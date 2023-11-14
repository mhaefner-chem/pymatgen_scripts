#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:21:39 2023

@author: bt308570
"""

# this program collects all results and compiles them into one file for each structure

import sys,os,importlib,glob,warnings,operator,math
from pymatgen.core import Structure,Composition
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


def energy(outcar):
    last_match = "free energy = 100000 eV"
    timing = "Total CPU : -100"

    
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

# get all configurations
configs = []
with open("supercell_coulomb_energy_l.txt","r") as file:
    for line in file:
        configs.append(line.split(".cif")[0])

calculations = [config for config in os.listdir(os.path.join(preopt_dir)) if os.path.isdir(os.path.join(preopt_dir,config))]
structures = {}


#specify the different calculation types
e_types = ["Q","fullAI","SP","fastAI","fastSP","SZP_SP","SZP_OPT","DZP_SP","DZP_OPT","ML"]#,"slowAI","slowSP"] # coulomb, sp, ML, partopt
e_min = [1] * len(e_types)
reference = "slowAI"

timings = {}
e_all = {}
run_success = {}
for config in configs:
    e_all[config] = [1] * len(e_types)
    timings[config] = [-1] * len(e_types)
    run_success[config] = [-1] * len(e_types)

# read in coulomb reference energies
with open("supercell_coulomb_energy_l.txt","r") as file:
    for line in file:
        key = line.split(".cif")[0]
        temp = line.split("\t")[1]
        index = e_types.index("Q")
        e_all[key][index] = float(temp.split(" ")[0])
        e_min[index] = min(e_all[key][index],e_min[index])


# get #atoms, composition, and Z from one of the structures (identical among all)
structure = Structure.from_file(config+"/POSCAR")
composition = Composition(structure.composition)
Z = composition.get_reduced_composition_and_factor()[1]
atoms = structure.num_sites

for calculation in calculations:
    config = calculation.split("_gpaw")[0]
    config = config.split("_fast")[0]
    config = config.split("_slow")[0]
    os.chdir(calculation)
    if os.path.isfile("done") and os.path.isfile("OUTCAR"):
        
        # get full ab-initio results
        if "fast" in calculation:
            index = e_types.index("fastAI")
        elif "slow" in calculation:
            index = e_types.index("slowAI")
        else:
            index = e_types.index("fullAI")
        
        final_energy,timing = energy("OUTCAR")
        if final_energy > 10000:
            print(calculation)
        run_success[config][index] = 0
        timings[config][index] = timing
        
        e_all[config][index] = final_energy/Z
        e_min[index] = min(e_all[config][index],e_min[index])
        
        # energy after single-point calculation
        
        if "fast" in calculation:
            index = e_types.index("fastSP")
        elif "slow" in calculation:
            index = e_types.index("slowSP")
        else:
            index = e_types.index("SP")
        e_all[config][index] = energy("SCFCONV/OUTCAR")[0]/Z
        run_success[config][index] = 0
        timings[config][index] = energy("SCFCONV/OUTCAR")[1]
        e_min[index] = min(e_all[config][index],e_min[index])
        
        
        # energy after at max 10 optimization steps
        if "partOpt" in e_types:
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
            run_success[config][index] = 0
            e_min[index] = min(e_all[config][index],e_min[index])
        
        # read in ML if it exists
        if os.path.isdir("ML") and "ML" in e_types:
            index = e_types.index("ML")
            e_all[config][index] = energy("ML/OUTCAR")[0]/Z
        e_min[index] = min(e_all[config][index],e_min[index])
        
    elif os.path.isfile("SCFCONV/OUTCAR"):
        if "fast" in calculation:
            index = e_types.index("fastAI")
        elif "slow" in calculation:
            index = e_types.index("slowAI")
        else:
            index = e_types.index("fullAI")
        run_success[config][index] = 1
        if "fast" in calculation:
            index = e_types.index("fastSP")
        elif "slow" in calculation:
            index = e_types.index("slowSP")
        else:
            index = e_types.index("SP")
        run_success[config][index] = 1
        
    elif os.path.isfile("SCFCONV/done") and os.path.isfile("SCFCONV/gpaw_sp.py"):
    # elif for GPAW results
    # read out results and timings in gpaw_e and gpaw_t
        if "dzp" in calculation:
            index = e_types.index("DZP_SP")
        elif "szp" in calculation:
            index = e_types.index("SZP_SP")
        
        i = 0
        with open("SCFCONV/GPAW.out", 'r') as file:
            for line in file:
                if "Extrapolated:" in line:
                    e = float(line.split(":")[1])/Z
                elif "Total:" in line:
                    timing = float(line.split()[1])
                    timings[config][index] = timing
                    run_success[config][index] = 0
        
        e_all[config][index] = e
        e_min[index] = min(e_all[config][index],e_min[index])
        
        if os.path.isfile("done") and os.path.isfile("gpaw_opt.py"):
        # elif for GPAW results
        # read out results and timings in gpaw_e and gpaw_t
            if "dzp" in calculation:
                index = e_types.index("DZP_OPT")
            elif "szp" in calculation:
                index = e_types.index("SZP_OPT")
            
            i = 0
            with open("GPAW.out", 'r') as file:
                for line in file:
                    if "Extrapolated:" in line:
                        e = float(line.split(":")[1])/Z
                    elif "Total:" in line:
                        run_success[config][index] = 0
                        timing = float(line.split()[1])
                        timings[config][index] = timing
                    if "TIME LIMIT" in line or "LineSearch failed" in line:
                        run_success[config][index] = 1
            
            e_all[config][index] = e
            e_min[index] = min(e_all[config][index],e_min[index])    
        elif os.path.isfile("gpaw_opt.py"):
            if "dzp" in calculation:
                index = e_types.index("DZP_OPT")
            elif "szp" in calculation:
                index = e_types.index("SZP_OPT")
            
            temp = glob.glob("*sum")
            for i in temp:
                with open(i, 'r') as file:
                    for line in file:
                        if "TIME LIMIT" in line or "LineSearch failed" in line:
                            run_success[config][index] = 1

    
    # else:
    #     e_all[config] = [1000.0] * len(e_types)
    #     if os.path.isfile("SCFCONV/done"):
    #         if os.path.isfile("SCFCONV/vasprun.xml"):
    #             xml = Vasprun("SCFCONV/vasprun.xml")
    #         else:
    #             print("Missing xml for",config)
    #         # if xml.converged_electronic == False:
    #         #     print("The SCF convergence failed!")                
    #     else:
    #         print(config,"Config not done!")
    
    
    os.chdir(os.path.join(preopt_dir))



# DO REGRESSION STUFF

# regression for coulomb prediction vs. optimized
# filter energies

e_proper = {}
for key,e in e_all.items():
    if e[0] < 0:
        e_proper[key] = e



x_lin = np.linspace(0, 2, 100)
ref = []
index = e_types.index(reference)

for key, values in e_proper.items():
    ref.append(values[index]/e_min[index])

fig, ax = plt.subplots()

size = 25
linewidth = 0.5    
r2_all = [-1] * len(e_types)

for i in range(len(e_types)):
    x = []
    y = []
    for key, values in e_proper.items():
        if values[i] < 0.5 and values[index]/e_min[index] > 0.97:
            y.append(values[index]/e_min[index])
            x.append(values[i]/e_min[i])#*(values[index]/e_min[index]))
    if x != []:
        slope, intercept, r, p, std_err = stats.linregress(x,y)
        print(e_types[i],slope,intercept,e_min[i],e_min[index])
        r2_all[i] = r**2
        y_lin = slope * x_lin + intercept
        plt.scatter(x,y,label=e_types[i]+", r² = "+"{:6.3f}".format(r**2),edgecolors="#000000",s=size,linewidth=linewidth)
        plt.plot(x_lin, y_lin)



ax.grid(zorder=0,linestyle="--",alpha=0.5)
plt.ylim(ymin=0.975, ymax=1.005)
plt.xlim(xmin=0.92, xmax=1.005)
plt.legend()
ax.set_xlabel('rel. Energy/Predictors')
ax.set_ylabel('rel. Energy/PBEsol')
 
# plt.show() 


# results

# timing parameters

max_timing = [0] * len(e_types)
min_timing = [0] * len(e_types)
avg_timing = [0] * len(e_types)


for i in range(len(e_types)):
    counter = 0
    sum = 0
    for config in timings:
        if timings[config][i] > 0:
            counter += 1
            sum += timings[config][i]
        if counter > 0:
            avg_timing[i] = sum/counter
        else:
            avg_timing[i] = -1
        
print("{:^10}".format("Method"),"{:^6}".format("r²"),"{:^8}".format("a.m. t/s"),"{:^8}".format("speed-up"),"{:^6}".format("r² > Q?"))

index_ref = e_types.index(reference)
index_Q = e_types.index("Q")
for i in range(len(e_types)):
    if r2_all[i] > 0:
        if avg_timing[i] > 0:
            ratio = avg_timing[index_ref]/avg_timing[i]
        else:
            ratio = -1
        if r2_all[i] < r2_all[index_Q]:
            better = "{:6}".format("no")
        else:
            better = "{:6.3f}".format(r2_all[i])
        print("{:10}".format(e_types[i]),"{:6.3f}".format(r2_all[i]),"{:8.0f}".format(avg_timing[i]),"{:8.1f}".format(ratio),better)

print("################################")

successful = [0] * len(e_types)
for i in range(1,len(e_types)):
    counter = [0] * 2
    for config in e_all:
        if run_success[config][i] == -1:
            counter[0] += 1

        elif run_success[config][i] == 1:
            counter[1] += 1
        successful[i] = len(e_all)-counter[0]-counter[1]
    print("{:10}".format(e_types[i]),":",str(counter[0])+"/"+str(len(e_all)),
          "missing,",str(counter[1])+"/"+str(len(e_all)),"failed")
    for config in e_all:
        if run_success[config][i] != 0 and successful[i] > 0:
            print("check",config,e_types[i])

print("################################")



with open("results.dat", 'w') as file:
    for i in range(len(e_types)):
        file.write("{:10}".format(e_types[i]))
    file.write("\n")
    for config in configs:
        for i in range(len(e_types)):
            file.write("{:10.3f}".format(e_all[config][i]))
        file.write("{:^10}".format(config.split("_")[1]))
        file.write("\n")
    file.write("timing_average")
    for i in range(len(e_types)):
        file.write("{:9.1f}".format(avg_timing[i]))
    file.write("\n")
    file.write("successful_calculations")
    for i in range(len(e_types)):
        file.write("{:5.0f}".format(successful[i]))
    file.write("\n")
    file.write("e_min:")
    for i in range(len(e_types)):
        file.write("{:10.3f}".format(e_min[i]))
    file.write("\n")
    
# sorted energy analysis

e_sorted = {}
for i in range(len(e_types)):
    tmp = []
    k = 0
    for config in configs:
        k += 1
        tmp.append([k,e_all[config][i]])
    e_sorted[i] = sorted(tmp, key=operator.itemgetter(1))

ratios = [0.01,0.05,0.1,0.25] #[0.01,0.05,0.10,0.20,0.25,0.5,1.0]

for ratio in ratios:
    scope = math.ceil(ratio*len(e_all))
    print("%%%%%%\n"+str(scope)+"\n%%%%%")
    ref = []
    e_truncated = {}
    overlap = [0] * len(e_types)
    
    for i in range(scope):
        ref.append(e_sorted[index][i][0])
        
    for i in range(len(e_types)):
        tmp = []
        for j in range(scope):
            tmp.append(e_sorted[i][j][0])
        e_truncated[i] = tmp
        
    for i in range(len(e_types)):
        for j in range(scope):
            if ref[j] in e_truncated[i]:
                overlap[i] += 1
        print("{:12}".format(e_types[i])," - overlap: ","{:6.2f}".format(overlap[i]/scope*100.0E0)," %")
            
    

sys.exit()

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
