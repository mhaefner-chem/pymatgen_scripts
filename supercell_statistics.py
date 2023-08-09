#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 14:18:04 2023

@author: bt308570
"""

# this program does statistics

import sys,os,importlib,glob,warnings
from pymatgen.core import Structure,Composition
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


# beginning of main program


sys.stdout.flush()
warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
warnings.filterwarnings("ignore", message="your POTCAR is corrupted")
print("Collecting Results.")

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
os.chdir(preopt_dir)

spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
directory_tools = importlib.util.module_from_spec(spec)
sys.modules["name"] = directory_tools
spec.loader.exec_module(directory_tools)        

data = glob.glob("*.dat")

for item in data:
    with open(item,"r") as file:
        for line in file:
            if "Q" in line:
                e_types = line.split()
                break


norm = []
timing_sets = []

for item in data:
    e_dat = []
    e_min = [1] * len(e_types)
    timings = [0.0] * len(e_types)
    successes = [0] * len(e_types)

    with open(item,"r") as file:
        for line in file:
            if "Q" in line:
                continue
            elif "timing_average" in line:
                for i in range(1,len(e_types)+1):
                    timings[i-1] = float(line.split()[i])
            elif "successful_calculations" in line:
                for i in range(1,len(e_types)+1):
                    successes[i-1] = int(line.split()[i])
            else:
                e = line.split()
                for i in range(len(e)):
                    e[i] = float(e[i])
                e_dat.append(e)
        timing_sets.append([timings,successes])
        
        for i in range(len(e_types)):
            for values in e_dat:
                e_min[i] = min(values[i],e_min[i])
        
        for values in e_dat:
            temp = []
            for i in range(len(e_types)):
                if values[i] < 0:
                    temp.append(values[i]/e_min[i])
                else:
                    temp.append(2)
            norm.append(temp)
                    
# timings

total_calculations = [0] * len(e_types)
total_average = [0.0] * len(e_types)

for set in timing_sets:
    for i in range(len(e_types)):
       total_calculations[i] += set[1][i]

for set in timing_sets:
    for i in range(len(e_types)):
        if total_calculations[i] != 0:
            total_average[i] += set[0][i] * set[1][i]/total_calculations[i]
        else:
            total_average[i] += -1
            
index = e_types.index("fullAI")
print("{:10}".format("type"),"{:9}".format("t_avg/s"))
for i in range(len(e_types)):
    if not "Q" in e_types[i]:
        print("{:10}".format(e_types[i]),"{:9.1f}".format(total_average[i]),"{:6.2f}".format(total_calculations[i]/len(norm)*100.0),"% success",
              "{:5.1f}".format(total_average[index]/total_average[i]),"x faster than fullAI")

# regressions and plotting
        
x_lin = np.linspace(0, 2, 100)
index = e_types.index("fullAI")
fig, ax = plt.subplots()

size = 10
linewidth = 0.25    
r2_all = [-1] * len(e_types)

# plot_list = ["Q","fullAI","SP","fastAI","fastSP","SZP_SP","SZP_OPT","DZP_SP","DZP_OPT"]
plot_list = ["Q","SP","fastSP","SZP_SP","DZP_SP","fullAI"] #,"fastAI","SZP_OPT","DZP_OPT","fullAI"] #,"fastSP","SZP_SP","DZP_SP","fullAI"]


e_label = [""] * len(e_types)

for i in range(len(e_types)):
    if e_types[i] == "fullAI":
        e_label[i] = "PBE/Acc-opt"
    elif e_types[i] == "Q":
        e_label[i] = "Ewald"
    elif e_types[i] == "SP":
        e_label[i] = "PBE/Acc-SP"
    elif e_types[i] == "fastSP":
        e_label[i] = "PBE/Norm-SP"
    elif e_types[i] == "SZP_SP":
        e_label[i] = "GPAW/szp-SP"
    elif e_types[i] == "DZP_SP":
        e_label[i] = "GPAW/dzp-SP"
    elif e_types[i] == "fastAI":
        e_label[i] = "PBE/norm-opt"
    elif e_types[i] == "SZP_OPT":
        e_label[i] = "GPAW/szp-opt"
    elif e_types[i] == "DZP_OPT":
        e_label[i] = "GPAW/dzp-opt"
    else:
        e_label[i] = e_types[i]
        

for item in plot_list:
    i = e_types.index(item)
    x = []
    y = []
    for values in norm:
        if values[i] < 1.5 and values[index] < 1.5 and values[index] > 0.98:
            y.append(values[index])
            x.append(values[i])
    
    if "Q" in item:
        speedup = "--"
    else:
        speedup = "{:5.1f}".format(total_average[index]/total_average[i])
    
    slope, intercept, r, p, std_err = stats.linregress(x,y)
    r2_all[i] = r**2
    y_lin = slope * x_lin + intercept
    plt.scatter(x,y,alpha=0.8,marker="o",
                label=e_label[i]+", r²:"+"{:6.3f}".format(r**2)+", t/"+r'$\bar{t}$'+" "+speedup,
                edgecolors="#000000",s=size,linewidth=linewidth)
    plt.plot(x_lin, y_lin)



ax.grid(zorder=0,linestyle="--",alpha=0.5)
plt.ylim(ymin=0.978, ymax=1.0025)
plt.xlim(xmin=0.935, xmax=1.0025)
plt.legend()
plt.title("Predictive performance for Al:Na$_2$ZnCl$_4$")
ax.set_xlabel('rel. Energy/Predictors')
ax.set_ylabel('rel. Energy/PBE-D3(BJ)/500 eV/Accurate')
# ax.set_aspect(2)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
 
plt.show() 
fig.savefig('supercell_statistics.png', format='png',bbox_inches="tight",dpi=600)

index = e_types.index("DZP_SP")               
counter = 0
for i in range(len(norm)):
    if norm[i][index] > 0.995:
        counter += 1
print(counter,len(norm))

sys.exit()
            
