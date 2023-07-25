#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 12:04:29 2023

@author: bt308570
"""



import sys,os,time
from glob import glob
import importlib.util
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
    
    
    max_jobs = 120
    
    sys.stdout.flush()
    
    file_dir = os.getcwd()
    
    
    print("Running script to carry out supercell calculations.")
    # locate head directory
    workdir = find_topdir()
    os.chdir(workdir)
    
    # import the directory tools 
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)   
    
    preopt_dir = workdir+"/SUPERCELL/PREOPT"
    outer_loop = 0
    while outer_loop < 200:
        outer_loop += 1
        print(outer_loop)
        os.chdir(preopt_dir)
        dirs = [dir for dir in os.listdir(preopt_dir) if os.path.isdir(os.path.join(preopt_dir, dir))]
        
    
        
        process = subprocess.Popen(['squeue',"-u","bt308570"],
                     stdin=subprocess.DEVNULL,
                     stdout=subprocess.PIPE,
                     start_new_session=False)
        out, err = process.communicate()
        jobs = out.decode().count("bt308570")
        
        checks = ["done","run"]
        break_loop = False
        
        for dir in dirs:
            # check if doped structure
            if "base" in dir:
                if break_loop == True:
                    break
                os.chdir(dir)
                a = ""
                for i in range(0,len(dir)):
                    a += "*"
                    
                print(a)
                print(dir)
                print(a)
                for subdir in glob("supercell*"):
                    if os.path.isdir(subdir):
                        os.chdir(subdir)
                        if not any(os.path.isfile("SCFCONV/"+item) for item in checks) and jobs < max_jobs:
                            print("Calculation set up",dir,subdir)
                            # vsub
                            name = glob("*.out")[0].split(".out")[0]
                            ref_file = glob("*cif")[0]
                            if "gpaw" in subdir:
                                print("GPAW found")
                                process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",ref_file,"-n",name,"-p","--gpaw"],
                                    stdin=subprocess.DEVNULL,
                                    stdout=open(name+'.out', 'w'),
                                    stderr=subprocess.STDOUT,
                                    start_new_session=True,
                                    preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                                time.sleep(1)
                                print("Job:",str(jobs)+"/"+str(max_jobs))
                                jobs += 1
                            else:
                                process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",ref_file,"-n",name,"-p"],
                                    stdin=subprocess.DEVNULL,
                                    stdout=open(name+'.out', 'w'),
                                    stderr=subprocess.STDOUT,
                                    start_new_session=True,
                                    preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                                time.sleep(1)
                                print("Job:",str(jobs)+"/"+str(max_jobs))
                                jobs += 1
                        elif jobs > max_jobs-1:
                            print("Job quota met, not starting any new jobs.")
                            break_loop = True
                            break
                        
            
                        os.chdir(os.path.join(preopt_dir,dir))
                os.chdir(preopt_dir)
        time.sleep(1200)
            
    

    

    