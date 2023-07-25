#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:26:40 2023

@author: bt308570
"""


import sys,os,time
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
    
print("RUNNING REPEATER!!!")
i = 0
while i < 40:
    
    i = i + 1
    print("Step ","{:03}".format(i))
    subprocess.Popen(['python3','-u', workdir+'/bin/env_run.py'],
                     stdin=subprocess.DEVNULL,
                     stdout=open('cycle.out', 'a'),
                     stderr=subprocess.STDOUT,
                     start_new_session=False,
                     preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    time.sleep(21600)