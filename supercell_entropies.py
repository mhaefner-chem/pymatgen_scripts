#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 10:51:26 2023

@author: bt308570
"""

import math,sys
import numpy as np

def calc_entropy(e,m,c,T):
    
    # compute the partition function
    pf = 0
    k_B = 8.617333E-5 # eV/K
    e_min = min(e)
    boltzmann = []
    calc = [True] * len(e)
    
    # precalculate exponential terms
    for i in range(len(e)):
        boltzmann.append(math.exp(-(e[i]-e_min)/(k_B*T)))
        if boltzmann[i] < 1E-12:
            calc[i] = False
    
    for i in range(len(e)):
        pf = pf + m[i] * boltzmann[i]
        
    # p_i = boltzmann[i]/pf
    S = 0
    for i in range(len(e)):
        if calc[i]:
            S = S - k_B * m[i] * boltzmann[i]/pf * math.log(boltzmann[i]/pf)
    
    for i in range(len(e)):
        if calc[i]:
            if -k_B * m[i] * boltzmann[i]/pf * math.log(boltzmann[i]/pf)*1000.0 > 1/T:
                print("Config {:>7s}: {:10.6f} meV/K at energy {:6.3f} eV".format(c[i],-k_B * m[i] * boltzmann[i]/pf * math.log(boltzmann[i]/pf)*1000.0,e[i]-e_min))
    return S

def read_configs(file):
    
    energies = []
    multiplicities = []
    configs = []
    
    with open(file, 'r') as f:
        for line in f:
            if not "ID" in line:
                configs.append(line.split()[0])
                energies.append(float(line.split()[1]))
                multiplicities.append(int(line.split()[2]))
            
    
    return energies, multiplicities, configs
        

# beginning of main program
if __name__ == '__main__':
    energies = []
    multiplicities = []
    configs = []
    T = 1000 # temperature in Kelvin
    

    n = len(sys.argv)
    for i in range(0,n):
        # specify the cells and supercell size from cif
        if sys.argv[i] == "-i":
            print("Reading in the file "+sys.argv[i+1]+".")
            energies, multiplicities, configs = read_configs(sys.argv[i+1])
    
    
    n = 21
    energies = np.linspace(0.0,1.0,n)
    multiplicities = [1] * n
    configs = ["a"] * n
    
    entropy = calc_entropy(energies, multiplicities,configs,T)
    print("    T: {:10.3f}   K".format(T))
    print("   ΔS: {:10.3f} meV/K".format(entropy*1000))
    print("e_min: {:10.3f}  eV".format(min(energies)))
    print("-ΔS*T: {:10.3f}  eV".format(-entropy*T))
    print("g_min: {:10.3f}  eV".format(min(energies)-entropy*T))