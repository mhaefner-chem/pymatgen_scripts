#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:26:40 2023

@author: bt308570
"""


import sys,os
import subprocess,signal

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

    # locate head directory
    workdir = find_topdir()
    os.chdir(workdir)
    

    subprocess.Popen(['python3','-u', workdir+'/bin/converge_cifs.py'],
                     stdin=subprocess.DEVNULL,
                     stdout=open('run1.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=False)

    subprocess.Popen(['python3','-u', workdir+'/bin/phonocalc.py'],
                     stdin=subprocess.DEVNULL,
                     stdout=open('run2.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=False)
    
    subprocess.Popen(['python3','-u', workdir+'/bin/phonoresult.py'],
                     stdin=subprocess.DEVNULL,
                     stdout=open('run3.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=False)

    subprocess.Popen(['python3','-u', workdir+'/bin/result_collector.py'],
                     stdin=subprocess.DEVNULL,
                     stdout=open('run4.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=False)
