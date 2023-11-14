#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 14:18:04 2023

@author: bt308570
"""

# this program does statistics

import sys,os,importlib,glob,warnings,math,operator
from pymatgen.core import Structure,Composition
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
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

data = glob.glob("*Z*.dat")

for item in data:
    with open(item,"r") as file:
        for line in file:
            if "Q" in line:
                e_types = line.split()
                break


norm = []
relative = []
full = []
timing_sets = []
e_dat = []

for item in data:
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
    temp_norm = []
    temp_rel = []
    temp_full = []
    for i in range(len(e_types)):
        if values[i] < 0:
            temp_full.append(values[i]-e_min[i]) # or e_min[i]
            temp_norm.append(values[i]/e_min[i])
            if values[e_types.index("fullAI")] < 0:
                temp_rel.append(values[i]-e_min[i]-values[e_types.index("fullAI")]+e_min[e_types.index("fullAI")])
            else:
                temp_rel.append(200)
        else:
            temp_full.append(20000)
            temp_norm.append(2)
            temp_rel.append(200)
    full.append(temp_full)
    norm.append(temp_norm)
    relative.append(temp_rel)
                    
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
 
#plt.show() 
fig.savefig('supercell_statistics.png', format='png',bbox_inches="tight",dpi=600)

# plot relative energies

# a = sorted(e_dat, key=operator.itemgetter(1))

# fig, ax = plt.subplots()
# cmap = mpl.cm.tab10

# i = 0
# for item in relative:
#     i += 1
#     j = 0
#     for value in item:
#         j += 1
#         plt.scatter(j,value,marker="o",alpha=0.8,
#                 label="",s=size,linewidth=linewidth,color=cmap(j/len(e_types)))
#         if i == 1:
#             plt.scatter(0,-20,marker="o",
#                     label=e_label[j-1],s=size,linewidth=linewidth,color=cmap(j/len(e_types)))
            
# plt.legend()
# plt.ylim(ymin=-1.5, ymax=1.5)
# plt.xlim(xmin=-5, xmax=10)
# ax.grid(zorder=0,linestyle="--",alpha=0.5)
# plt.show()

# sorted energies analysis

e_sorted = {}
for i in range(len(e_types)):
    tmp = []
    k = 0
    for j in full:
        k += 1
        tmp.append([k,j[i]])
    e_sorted[i] = sorted(tmp, key=operator.itemgetter(1))

# ratios = [0.01,0.05,0.10,0.20,0.25,0.5,1.0]
ratios = []
#ratios = np.linspace(0.01,1.0,num=math.ceil(len(full)/5))

overlap_array = {}

# total overlap at all points
for ratio in ratios:
    scope = math.ceil(ratio*len(full))
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
    
    tmp = []
    for i in range(len(e_types)):
        
        for j in range(scope):
            if ref[j] in e_truncated[i]:
                overlap[i] += 1
        print("{:12}".format(e_types[i])," - overlap: ",overlap[i]/scope*100.0E0," %")
        tmp.append(overlap[i]/scope)
    overlap_array[ratio] = tmp
    
# write results to file
with open("overlaps.dat", 'w') as file:
    for i in range(len(e_types)):
        file.write(e_types[i]+" ")
    file.write("\n")

    for ratio in ratios:
        scope = math.ceil(ratio*len(full))
        file.write("{:4}".format(scope))
        for i in range(len(e_types)):
            file.write("{:7.4f}".format(overlap_array[ratio][i]))
        file.write("\n")

# overlap for the N lowest structures
        
N = 25
overlap_array = {}
ratios = np.linspace(0.02,0.3,num=math.ceil(len(full)*0.23))

for i in range(len(full)):
    print(i+1,full[i][index])
for i in range(N):
    print(e_sorted[index][i])
print("Q")
for i in range(N):
    print(e_sorted[e_types.index("Q")][i])
print("SZP_SP")
for i in range(N):
    print(e_sorted[e_types.index("SZP_SP")][i])
print("SZP_OPT")
for i in range(N):
    print(e_sorted[e_types.index("SZP_OPT")][i])

#for i in range(N):
#    print(full[e_sorted[index][i][0]][index])
for ratio in ratios:
    scope = math.ceil(ratio*len(full))
    ref = []
    e_truncated = {}
    overlap = [0] * len(e_types)
    
    for i in range(N):
        ref.append(e_sorted[index][i][0])
    
    for i in range(len(e_types)):
        tmp = []
        for j in range(scope):
            tmp.append(e_sorted[i][j][0])
        e_truncated[i] = tmp
    
    tmp = []
    for i in range(len(e_types)):
        
        for j in range(scope):
            if e_truncated[i][j] in ref:
                overlap[i] += 1
        tmp.append(overlap[i]/N)
    overlap_array[ratio] = tmp
        
with open("lowest_overlaps.dat", 'w') as file:
    for i in range(len(e_types)):
        file.write(e_types[i]+" ")
    file.write("\n")

    for ratio in ratios:
        scope = math.ceil(ratio*len(full))
        file.write("{:4}".format(scope))
        for i in range(len(e_types)):
            file.write("{:7.4f}".format(overlap_array[ratio][i]))
        file.write("\n")



index = e_types.index("DZP_SP")               
counter = 0
for i in range(len(norm)):
    if norm[i][index] > 0.995:
        counter += 1
print(counter,len(norm))

sys.exit()
            
