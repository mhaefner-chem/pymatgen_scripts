#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 14:57:48 2023

@author: bt308570
"""

# make asymmetric
# get time estimate from first step
# give resolution in AA

from pymatgen.core import Structure,Species
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer,SpacegroupOperations
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.vasp import Poscar
import matplotlib.pyplot as plt
import numpy as np
import sys,time


def line(step,resolution):
    line = np.linspace(0.01,0.99,resolution)
    return line[step]

def distance(a,b):
    dist = ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
    return dist

n = len(sys.argv)
for i in range(0,n):
    # specify the cell from cif, POSCAR, or CONTCAR
    if sys.argv[i] == "-c":
        print("Reading in the file "+sys.argv[i+1]+".")
        struc = Structure.from_file(sys.argv[i+1])
    elif sys.argv[i] == "-s":
        print("Substituted element: "+sys.argv[i+1])
        switch = int(sys.argv[i+1])



resolution = 0.25 # in AA
cutoff = 3 # eV
limit_low_cutoff = 0.85


res = [1,1,1]
real_res = [resolution,resolution,resolution]
steps = 1
for i in range(0,3):
    res[i] = int(np.ceil(struc.lattice.abc[i]/resolution))
    real_res[i] = struc.lattice.abc[i]/res[i]
    steps = res[i]*steps

print("Real resolution: ",real_res)

for site in struc.sites:
    if "Na" in site.species:
        site.species = "Na1+"
    elif "Zn" in site.species:
        site.species = "Zn2+"
    elif "Al" in site.species:
        site.species = "Al3++"
    elif "Cl" in site.species:
        site.species = "Cl1-"
    else:
        site.species = "H0"

combined_radii = {}
sub_ion = struc.sites[switch-1]
for site in struc.sites:
    combined_radii[site.specie.symbol] = site.specie.ionic_radius+sub_ion.specie.ionic_radius


time_large_loop = {}


e_tot = np.zeros((res))
for x in range(0,res[0]):
    time_large_loop[x] = time.thread_time()
    if x > 0:
        time_per_step = time_large_loop[x] - time_large_loop[x-1]
        print("Time for step/sec:  ","{:8.3f}".format(time_per_step))
        print("Remaining time/sec: ","{:8.3f}".format(time_per_step*(res[0]-x)))
    for y in range(0,res[1]):
        for z in range(0,res[2]):
            new_coords = [line(x,res[0]),line(y,res[1]),line(z,res[2])]
            struc.replace(switch-1,species="Na1+",coords=new_coords)
            neighbors = struc.get_neighbors(struc.sites[switch-1],r=3.5)
            for neighbor in neighbors:
                if neighbor.nn_distance < limit_low_cutoff*combined_radii[neighbor.specie.symbol]:
                    substitute = False
                    break
                else:
                    substitute = True
            if substitute == True:

                ES = EwaldSummation(struc)
                e_tot[x,y,z] = ES.total_energy
          
print("Total time/sec:     ","{:8.3f}".format(time_large_loop[res[0]-1]-time_large_loop[0]))

min_index = np.unravel_index(np.argmin(e_tot), e_tot.shape)
min_val = e_tot[min_index]


indices = np.where(e_tot <= min_val+cutoff)


stable_options = []


indices_3d = np.transpose(indices)

for i in indices_3d:
    k = [0,0,0]
    for j in range(0,3):
        k[j] = i[j]/res[j]    
    if e_tot[i[0],i[1],i[2]] < -1:
        stable_options.append([k,e_tot[i[0],i[1],i[2]]])
        struc.append("H0",k)
struc.remove_sites([switch-1])
poscar = Poscar(struc)
poscar.write_file("POSITIONS.vasp")
    
          
sys.exit()

            # coords = struc.lattice.get_cartesian_coords(fractional_coords=[line(x,res[0]),line(y,res[1]),line(z,res[2])])
                        

            
            # if len(struc.get_neighbors_in_shell(coords,r=0.0,dr=limit_low)) > 0:
            #     attrac_all[x,y,z] = 10000.0
            #     continue
            # else:
            #     neighbors = struc.get_neighbors_in_shell(coords,r=limit_low,dr=limit_high)
            #     for neighbor in neighbors:
            #         # print(neighbor.nn_distance,neighbor.specie.symbol,neighbor.specie.oxi_state)
            #         if neighbor.nn_distance < limit_low_cutoff*combined_radii[neighbor.specie.symbol]:
            #             local_attraction = local_attraction + 10000.0
            #         else:
            #             E_C = 1.43965E1*neighbor.specie.oxi_state*sub_species.oxi_state/neighbor.nn_distance
            #             local_attraction = local_attraction + E_C
                

                
            #     # local_attraction = len(struc.get_neighbors_in_shell(coords,r=limit_low,dr=limit_high))
            
            # attrac_all[x,y,z] = local_attraction
            



    
k_cart = [0.0,0.0,0.0]
l_cart = [0.0,0.0,0.0]
    
for i in range(0,len(stable_options)):
    for j in range(i+1,len(stable_options)):
        for k in stable_options[i][0]:
            for n in range(3):
                k_cart[n] = k[n]*struc.lattice.abc[n]
                
            for l in stable_options[j][0]:
                for n in range(3):
                    l_cart[n] = l[n]*struc.lattice.abc[n]
                    
                if distance(k_cart,l_cart) < sub_species.ionic_radius:
                    if stable_options[i][1] < stable_options[j][1]:
                        stable_options[j][1] = 20000.0
                    else:
                        stable_options[i][1] = 20000.0
                
                # print(k,l,distance(k_cart,l_cart))

print("")
print("Obtained vacant positions:")
pos = [1.0,1.0,1.0]
for i in range(0,len(stable_options)):
    if stable_options[i][1] < stab_limit:
        print("Position (frac):",end=' ')
        
        
        mini = 1.0
        index = -1
        
        for j in range(len(stable_options[i][0])):
            if stable_options[i][0][j][0] < mini:
                mini = stable_options[i][0][j][0]
                index = j
            
        for j in range(3):
            pos[j] = "{:5.3f}".format(stable_options[i][0][index][j])
            
            
        print(pos[0]+", "+pos[1]+", "+pos[2],end=',')
        testpos.append([stable_options[i][0][index],stable_options[i][1]])
        print(" local stability/eV:","{:8.3f}".format(stable_options[i][1])," , multiplicity: ",mult)
        for j in range(len(stable_options[i][0])):
            print("(",end="")
            for k in range(3):
                print("{:5.3f}".format(stable_options[i][0][j][k]),end=" ")
            print(")",end=" ")
            print("")

sys.exit()

combined_radii = {}

space_group = SpacegroupAnalyzer(struc)
print("Multiplicity:",len(space_group.get_space_group_operations()))
mult = len(space_group.get_space_group_operations())

for site in struc.sites:
    combined_radii[site.specie.symbol] =                    site.specie.ionic_radius+sub_species.ionic_radius





# neighborcount = np.zeros((resolution,resolution,resolution))
# neighborline = np.zeros((resolution,resolution))
# neighborplane = np.zeros((resolution))



for limit_high_cutoff in cutoffs:

    print("***********************************")
    print("Cutoff: ",limit_high_cutoff)
    print("***********************************")
    limit_low = combined_radii[min(combined_radii)]*limit_low_cutoff
    limit_high = combined_radii[max(combined_radii)]*limit_high_cutoff


    time_large_loop = {}


    attrac_all = np.zeros((res))
    for x in range(0,res[0]):
        time_large_loop[x] = time.thread_time()
        if x == 1:
            time_per_step = time_large_loop[x] - time_large_loop[x-1]
            print("Time for step/sec:  ","{:8.3f}".format(time_per_step))
            print("Remaining time/sec: ","{:8.3f}".format(time_per_step*(res[0]-x)))
            # print("Step:",x+1," of ",res[0])
        for y in range(0,res[1]):
            for z in range(0,res[2]):
                
                coords = struc.lattice.get_cartesian_coords(fractional_coords=[line(x,res[0]),line(y,res[1]),line(z,res[2])])
                            
                local_attraction = 0.0
                
                if len(struc.get_neighbors_in_shell(coords,r=0.0,dr=limit_low)) > 0:
                    attrac_all[x,y,z] = 10000.0
                    continue
                else:
                    neighbors = struc.get_neighbors_in_shell(coords,r=limit_low,dr=limit_high)
                    for neighbor in neighbors:
                        # print(neighbor.nn_distance,neighbor.specie.symbol,neighbor.specie.oxi_state)
                        if neighbor.nn_distance < limit_low_cutoff*combined_radii[neighbor.specie.symbol]:
                            local_attraction = local_attraction + 10000.0
                        else:
                            E_C = 1.43965E1*neighbor.specie.oxi_state*sub_species.oxi_state/neighbor.nn_distance
                            local_attraction = local_attraction + E_C
                    
    
                    
                    # local_attraction = len(struc.get_neighbors_in_shell(coords,r=limit_low,dr=limit_high))
                
                attrac_all[x,y,z] = local_attraction
                
    
    
    print("Total time/sec:     ","{:8.3f}".format(time_large_loop[res[0]-1]-time_large_loop[0]))
    
    min_index = np.unravel_index(np.argmin(attrac_all), attrac_all.shape)
    min_val = attrac_all[min_index]
    
    
    indices = np.where(attrac_all <= stab_limit)
    
    
    stable_options = []
    
    k = [0,0,0]
    indices_3d = np.transpose(indices)
    
    for i in indices_3d:
        # print(i/resolution,attrac_all[i[0],i[1],i[2]])
        points = []
        for j in range(0,3):
            k[j] = i[j]/res[j]
        for ops in space_group.get_symmetry_operations():
            newpoint = ops.operate(k)
            newpoint[newpoint < 0.0] += 1
            newpoint[newpoint > 1.0] += -1
            points.append(newpoint)
            #print(newpoint)
        
        stable_options.append([points,attrac_all[i[0],i[1],i[2]]])
        
    k_cart = [0.0,0.0,0.0]
    l_cart = [0.0,0.0,0.0]
        
    for i in range(0,len(stable_options)):
        for j in range(i+1,len(stable_options)):
            for k in stable_options[i][0]:
                for n in range(3):
                    k_cart[n] = k[n]*struc.lattice.abc[n]
                    
                for l in stable_options[j][0]:
                    for n in range(3):
                        l_cart[n] = l[n]*struc.lattice.abc[n]
                        
                    if distance(k_cart,l_cart) < sub_species.ionic_radius:
                        if stable_options[i][1] < stable_options[j][1]:
                            stable_options[j][1] = 20000.0
                        else:
                            stable_options[i][1] = 20000.0
                    
                    # print(k,l,distance(k_cart,l_cart))
    
    print("")
    print("Obtained vacant positions:")
    pos = [1.0,1.0,1.0]
    for i in range(0,len(stable_options)):
        if stable_options[i][1] < stab_limit:
            print("Position (frac):",end=' ')
            
            
            mini = 1.0
            index = -1
            
            for j in range(len(stable_options[i][0])):
                if stable_options[i][0][j][0] < mini:
                    mini = stable_options[i][0][j][0]
                    index = j
                
            for j in range(3):
                pos[j] = "{:5.3f}".format(stable_options[i][0][index][j])
                
                
            print(pos[0]+", "+pos[1]+", "+pos[2],end=',')
            testpos.append([stable_options[i][0][index],stable_options[i][1]])
            print(" local stability/eV:","{:8.3f}".format(stable_options[i][1])," , multiplicity: ",mult)
            for j in range(len(stable_options[i][0])):
                print("(",end="")
                for k in range(3):
                    print("{:5.3f}".format(stable_options[i][0][j][k]),end=" ")
                print(")",end=" ")
                print("")




for i in range(0,len(testpos)):
    for j in range(i+1,len(testpos)):
        for n in range(3):
            k_cart[n] = testpos[i][0][n]*struc.lattice.abc[n]
                
            for n in range(3):
                l_cart[n] = testpos[j][0][n]*struc.lattice.abc[n]
                    
                if distance(k_cart,l_cart) < sub_species.ionic_radius:
                    if testpos[i][1] < testpos[j][1]:
                        testpos[j][1] = 20000.0
                    else:
                        testpos[i][1] = 20000.0


print("***********************************")
print("Final suggestions:")
for i in range(0,len(testpos)):
    if testpos[i][1] < 10000.0:
        for j in range(3):
            pos[j] = "{:5.3f}".format(testpos[i][0][j])
            
        print(pos[0]+", "+pos[1]+", "+pos[2],end='\n')

sys.exit()
            
                        
# min_index = np.unravel_index(np.argmin(neighborcount), neighborcount.shape)
# max_index = np.unravel_index(np.argmax(neighborcount), neighborcount.shape)

# min_val = neighborcount[min_index]
# max_val = neighborcount[max_index]

# print("Minimum value:", min_val)
# print("Index of minimum value:", min_index)
# print("Maximum value:", max_val)
# print("Index of maximum value:", max_index)         
# print("Cutoff value:",min_val+0.1*(max_val-min_val))

# indices = np.where(neighborcount <= (min_val+0.05*(max_val-min_val)))
# indices_3d = np.transpose(indices)
# print("3D indices of values within 10% of the minimum value:")
# print(indices_3d/resolution)

# for z in range(0,resolution):
#     fig, ax = plt.subplots()
#     im = ax.imshow(neighborcount[:,:,z])
    
#     # Show all ticks and label them with the respective list entries
#     #ax.set_xticks()
#     # ax.set_yticks(np.arange(len(vegetables)), labels=vegetables)
    
#     # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#              rotation_mode="anchor")
    
#     # Loop over data dimensions and create text annotations.
#     # for i in range(0,resolution):
#     #     for j in range(0,resolution):
#     #          text = ax.text(j, i, "{:6.2f}".format(neighborcount[i,j,1]),
#                             #ha="center", va="center", color="w")
    
#     # ax.set_title("Harvest of local farmers (in tons/year)")
#     fig.tight_layout()
#     plt.show()
#     time.sleep(1)
