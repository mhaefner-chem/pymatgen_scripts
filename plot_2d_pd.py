# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 15:38:33 2023

@author: bt308570
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import sys,os,importlib,glob,math


def linear_function(a,b):
    m = (a[1]-b[1])/(a[0]-b[0])
    b = -(m*a[0] + a[1])
    return m,b

# settings

tmax = 1200.0
max_steps = 10

sys.stdout.flush()
print("Collecting Results.")

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
        
spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
directory_tools = importlib.util.module_from_spec(spec)
sys.modules["name"] = directory_tools
spec.loader.exec_module(directory_tools)        

 # move through the results
files = glob.glob(workdir+"/RESULTS/*_delta.txt")


pd_struc_name = []
pd_struc_deltaG = {}
pd_struc_deltaE = {}
pd_ref_states = {}
all_ref_states = []

# read in information about the various states
for file in files:
    thermo_g = []
    ref_states = {}

    with open(file,'r') as f:
        lines = f.readlines()
        name = lines[0]
        for i in range(1,len(lines)):
            if "T/K" in lines[i]:
                for j in range(i+1,len(lines)):
                    thermo_g.append(lines[j].split())
            elif "Delta_E_0" in lines[i]:
                pd_struc_deltaE[name] = lines[i].split()[1]
            elif "SG" in lines[i]:
                # determine the reference states for each structure
                ref_states[lines[i].split()[0]] = float(lines[i].split()[1])
                if not lines[i].split()[0] in all_ref_states:
                    all_ref_states.append(lines[i].split()[0])
                
    pd_struc_deltaG[name] = thermo_g
    pd_struc_name.append(name)
    pd_ref_states[name] = ref_states
            


# Define the x and y coordinates for the points

max_y = 0.0
min_y = 0.0

# do the binary PD plot for all binary structures
for i in range(0,len(all_ref_states)):
    for j in range(i+1,len(all_ref_states)):
        print(all_ref_states[i],all_ref_states[j])
        ref0 = [0.00, 0.00]
        ref1 = [100.00, 0.00]
        
        # Create a new figure
        fig, ax = plt.subplots()

        # temperature-based results with colormap

        cmap = mpl.cm.brg
        norm = mpl.colors.Normalize(vmin=0, vmax=tmax)
        ax.set_xlim([0, tmax])

        plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='horizontal', label='Temperature')




        # Connect the points with lines
        ax.plot([ref0[0], ref1[0]], [ref0[1], ref1[1]], '-', color='#000000', linewidth=1.5,zorder=0)
        ax.plot([ref0[0], ref1[0]], [ref0[1], ref1[1]], '-', color='#FFFFFF', linewidth=0.75,zorder=1)
        
        t_all = []
        p_all = []
        label = 0
        
        # add the results for the various structures
        for struc in pd_struc_name:
            total = 0
            for ref in pd_ref_states[struc]:
                total = total + pd_ref_states[struc][ref]

            if pd_ref_states[struc][all_ref_states[i]] > 0.1 and pd_ref_states[struc][all_ref_states[j]] > 0.1:
                
                divider = math.ceil(len(pd_struc_deltaG[struc])/max_steps)
                a = 0
                b = 0
                ratio = pd_ref_states[struc][ref]/total
                for thermo in pd_struc_deltaG[struc]:
                    a = a + 1
                    if float(thermo[0]) > tmax:
                        break
                    if not a%divider == 1:
                        continue
                    else:
                        t_all.append(thermo[0])
                    
                    
                    p_all.append([ratio*100.0,float(thermo[1]),thermo[0]])
                    label = label + 1
                    
                    ax.plot([ratio*100.0, ref1[0]], [float(thermo[1]), ref1[1]], '-', c="#000000", linewidth = 0.25,zorder=j-50)
                    ax.plot([ref0[0], ratio*100.0], [ref0[0],float(thermo[1])], '-', c="#000000", linewidth = 0.25,zorder=j-50)
                    
                    
                    
                    
                    
                    ax.scatter(ratio*100.0,float(thermo[1]), s=25, alpha = 0.333,label="", c="#FFFFFF",marker='o',edgecolors=cmap(float(thermo[0])/tmax),linewidth=1.0,zorder=j+50)
                    
                    if float(thermo[1]) > max_y:
                        max_y = float(thermo[1])
                    if float(thermo[1]) < min_y:
                        min_y = float(thermo[1]) 
                    
        for t in t_all:
            m_sub = []
            p_sub = []
            for p in p_all:
                if p[2] == t:
                    p_sub.append([p[0],p[1]])
            for p in p_sub:
                if p[1] < 0.0:
                    print(1)
            
            
        
        # Set the x and y axis limits
        ax.set_xlim([0, 100])
        ax.set_ylim([min_y*1.2, max_y*1.2])
        
        # Add labels for the x and y axes
        ax.grid(zorder=0,linestyle="--",alpha=0.5)
        ax.set_xlabel('Ratio '+all_ref_states[j].split("_")[0]+'[%]')
        ax.set_ylabel('Energy/atom [eV]')

        ax.legend()
        plt.title("Phase Diagram "+all_ref_states[i].split("_")[0]+"-"+all_ref_states[j].split("_")[0],fontsize=16,loc="center",pad=24.0)
        break
        
        
fig.savefig("PD_"+all_ref_states[i].split("_")[0]+"_"+all_ref_states[j].split("_")[0]+".svg", format='svg',bbox_inches="tight")
