#!/home/70/bt308570/anaconda3/envs/pymatgen/bin/python
from pymatgen.core import structure
from pymatgen.io.vasp.sets import MPScanRelaxSet, MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun,Outcar
import pymatgen.io.vasp.inputs as vasp_i
from pymatgen.io.vasp.inputs import Incar,Kpoints
import time,warnings,math

import sys,os,subprocess,shutil,importlib

# sets up a basic INCAR for an R2SCAN-D4 optimization

# detects the head directory
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

# creates the VASP inputs (INCAR, POTCAR, POSCAR)
def input_maker(param = 0,max_steps = 99):
    if param == 0:
        vasp_input = MPScanRelaxSet(struc)
        vasp_input.user_incar_settings = settings.incar_settings("SCANopt")
        vasp_input.user_incar_settings["NSW"] = max_steps
    elif param == 1:
        vasp_input = MPRelaxSet(struc,user_potcar_functional="PBE_54")
        vasp_input.user_incar_settings = settings.incar_settings("PSopt")
        vasp_input.user_incar_settings["NSW"] = max_steps * 2
    accuracy = 12.0
    vasp_input.user_incar_settings["NGX"] = math.ceil(struc.lattice.abc[0] * accuracy)
    vasp_input.user_incar_settings["NGY"] = math.ceil(struc.lattice.abc[1] * accuracy)
    vasp_input.user_incar_settings["NGZ"] = math.ceil(struc.lattice.abc[2] * accuracy)
        
    vasp_input.write_input(".")
    incar = Incar.from_file("INCAR")
    if "GGA" in incar and "METAGGA" in incar:
        incar.pop("METAGGA")
    incar.write_file("INCAR")
    return

# runs a VASP calculation based on a given vsub and directory
def vsub(name,dir):
        
    if os.path.isfile("done") == True:
        print("Calculation already done! Skipping to analysis.")
        #result = subprocess.run(["ls -alrt"], shell=True, stdout=subprocess.PIPE)
        #print(result.stdout.decode())
        result = Vasprun("vasprun.xml")
        return result.final_energy
    elif os.path.isfile("run") == True:
        print("Calculation is still running! Stopping script.")
        sys.exit()
    else:
        print("Calculation is set up.")
        workdir = os.getcwd()
        os.chdir(dir)
        if os.path.exists("POSCAR"):
            shutil.copyfile("POSCAR",name+".inp")
            
        # do an additional loop to catch the ZBRENT exception
        i = 0
        while i == 0:
            shutil.copyfile("INCAR","run")  
            print(subprocess.run(["/home/70/bt308570/bin/vsub_py "+name+".inp"], shell=True, stdout=subprocess.PIPE)) 
            sleeper()
            
            j = 1
            with open('OUTCAR','r') as file:
                for line in file:
                    if 'ZBRENT: fatal error in bracketing' in line:
                        shutil.copyfile("CONTCAR",name+".inp")
                        os.remove("done")
                        print("ZBRENT error encountered!")
                        j = 0
            i = j
        result = Vasprun("vasprun.xml")
        os.chdir(workdir)
        return result.final_energy 

def gsub(gpaw_mode,name):
    if os.path.isfile("done") == True:
        print("Calculation already done! Skipping to analysis.")
    elif os.path.isfile("run") == True:
        print("Calculation is still running! Stopping script.")
        sys.exit()
    else:
        print("Calculation is set up.")
        shutil.copyfile(name+".inp","run")  
        shutil.copyfile(name+".inp","POSCAR") 
        print(subprocess.run(["/home/70/bt308570/bin/gsub "+gpaw_mode+".py"], shell=True, stdout=subprocess.PIPE)) 
        sleeper()
        return

# sets up a sleep loop for the script to wait for the calculation to finish, runs for about a week before giving up
def sleeper():
    for i in range(1,10000):
        # print("This job has been running for: ",i," min")
        time.sleep(60)
        if os.path.isfile("done") == True:
            return
    print("Couldn't wait until job was finished. Giving up.")
    sys.exit()

# help function to output the commands and everything
def program_help():
    print("**********\nThis script can set up simple VASP optimizations with R2SCAN-D4.\nIn its most basic version, it only requires a .cif-file, POSCAR, CONTCAR, or .vasp-file as initial geometry\n    optimizer.py -c POSCAR\n")
    print("A custom INCAR can be supplied with '-i INCAR' and a custom name with '-n NAME'.")
    print("This help dialogue is displayed with either '-h' or '--help'.")
    print("A full optimization and convergence of ENCUT and KPOINTS can be achieved via '--fullopt'.")
    print("In case you only want to generate an INCAR without any calculations, e.g., to specify MAGMOM, you can use '--ionly'")
    print("**********")
    sys.exit()

# converge the k-points
def convergence(name,values,prop):
    workdir = os.getcwd()
    energy = {}
    xml = {}
    for i in values:
        tmpdir = workdir+"/"+name+"_"+prop+"_"+str(i)
        input_copy(name+"_"+prop+"_"+str(i),tmpdir)
        os.chdir(tmpdir)
        incar = vasp_i.Incar.from_file("INCAR")
        incar[prop] = i
        incar["NSW"] = 0
        incar.write_file("INCAR")
        
        print("Performing single-point calculation for "+prop+" = "+str(i)+".")
        energy[i] = vsub(name+"_"+prop+"_"+str(i),tmpdir)
        xml[i] = Vasprun("vasprun.xml")
        os.chdir(workdir)
    return energy,xml

# copies an existing input to a new directory
def input_copy(name,dir):
        if not os.path.exists(dir):
            os.makedirs(dir)
        if os.path.exists("CONTCAR"):
            shutil.copyfile("CONTCAR",dir+"/"+name+".inp")
        else:
            shutil.copyfile("POSCAR",dir+"/"+name+".inp")
        shutil.copyfile("INCAR",dir+"/INCAR")
        shutil.copyfile("POTCAR",dir+"/POTCAR")

# beginning of main program
if __name__ == '__main__':
    sys.stdout.flush()
    param = 0
    preopt = 0
    gpaw  = False
    warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
    n = len(sys.argv)
    for i in range(0,n):
        # specify the cell from cif, POSCAR, or CONTCAR
        if sys.argv[i] == "-c":
            print("Reading in the file "+sys.argv[i+1]+".")
            struc = structure.Structure.from_file(sys.argv[i+1])
        # should an external INCAR be considered?
        elif sys.argv[i] == "-i":
            print("Reading in the file "+sys.argv[i+1]+" as INCAR.")
            incar = vasp_i.Incar.from_file(sys.argv[i+1])
        # full optimization of geometry, k-points, and cutoff
        elif sys.argv[i] == "--fullopt":
            print("Convergence of k-points requested.")
            param = 2            
        elif sys.argv[i] == "--ionly":
            param = 3
        elif sys.argv[i] == "--gpaw":
            gpaw = True
        elif sys.argv[i] == "-n" or sys.argv[i] == "--name":
            name = sys.argv[i+1]
        elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
            program_help()
        elif sys.argv[i] == "-p" or sys.argv[i] == "--preopt":
            preopt = 1

if not "struc" in globals():
    print("No geometry data supplied with '-c'. Stopping script.")
    sys.exit()
if not "name" in globals():
    print("No name supplied with '-n'. Stopping script.")
    sys.exit()
    
# imports the settings.py file for the input generation
spec = importlib.util.spec_from_file_location("settings",find_topdir() + "/bin/settings.py")
settings = importlib.util.module_from_spec(spec)
sys.modules["name"] = settings
spec.loader.exec_module(settings)

# creates or retains the INCAR, preopt switches between settings detailed in "settings.py"
if "incar" in globals():
    input_maker(preopt)
    incar.write_file('INCAR')
else:
    input_maker(preopt)

# program stops if only the input was requested, but no actual calculation
if param == 3:
    sys.exit()

# get working directory
workdir = os.getcwd()

# preliminary SCF convergence test, also gives a conservative estimate on the maximum amount of optimization cycles possible in 24 hrs
scfdir = workdir+"/SCFCONV"
input_copy(name,scfdir)
if gpaw == True:
    shutil.copy("KPOINTS", "SCFCONV/KPOINTS")
os.chdir(scfdir)
if gpaw == False:
    incar = vasp_i.Incar.from_file("INCAR")
    incar["NSW"] = 0
    incar.write_file("INCAR")
    print("Running preliminary SCF convergence test.")
    vsub(name,scfdir)
    result = Vasprun("vasprun.xml")
    outcar = Outcar("OUTCAR")
    max_steps = 80000.0/outcar.run_stats["Elapsed time (sec)"]
else:
    if os.path.isfile("KPOINTS"):
        kpoints = Kpoints.from_file("KPOINTS")
        new_kpts = "("
        j = 0
        for i in kpoints.kpts[0]:
            new_kpts += str(i)
            j += 1
            if j == len(kpoints.kpts[0]):
                new_kpts += ")"
            else:
                new_kpts += ","
        
    with open(find_topdir()+"/bin/gpaw_sp.py","rt") as fin:
        with open("gpaw_sp.py", "wt") as fout:
            for line in fin:
                if "kpts={" in line:
                    fout.write(line.replace("(4,4,4)",new_kpts))
                elif "basis=" in line:
                    fout.write(line.replace("dzp","szp"))
                # if "txt='GPAW.out'" in line:
                #     fout.write(line.replace("GPAW",name))
                else:
                    fout.write(line)
    
    
    print("Running preliminary SCF convergence test with GPAW.")
    gpaw_mode = "gpaw_sp"
    gsub(gpaw_mode,name)
    sys.exit()


# kills job if the convergence failed
if result.converged_electronic == False:
    print("The initial SCF convergence failed! Please check the validity of the input files and provide an alternative INCAR if necessary.")
    sys.exit()
    
# runs vsub_py for optimization  
os.chdir(workdir)
input_maker(preopt,int(max_steps))

a = 0
while True:
    energy = vsub(name,workdir)
    result = Vasprun("vasprun.xml")
    a = a + 1
    # repeats optimization if not yet converged
    if result.converged_ionic == False:
        print("Optimization not converged on step",a,". Updating POSCAR to do another optimization.")
        os.remove("done")
        shutil.copyfile("POSCAR","POSCAR_"+str(a))
        shutil.copyfile("CONTCAR","POSCAR")
    elif result.converged_ionic == True:
        break
    elif result.converged_ionic == False and a > 10:
        print("Optimization failed even after 10 restarts. Please check system manually.")
        sys.exit()

# converges k-points if requested
if param == 2:
    # # encut first
    # energies_encut = convergence(name,[300,400,500,600,700,800],"ENCUT")
    # for i in energies_encut:
    #     if abs(energies_encut[i]-energies_encut[800]) < 0.01:
    #         print("This structure is converged at a cutoff energy of",i,"eV.")
    #         conv_cut = i
    #         break
       # else:
       #     print("ENCUT for this structure could not be converged with the given set of energies. Please check and do some additional testing.")
       
    energies_kpoint, xmls = convergence(name,[0.5,0.4,0.3,0.2],"KSPACING")
    for i in energies_kpoint:
        if i == 0.2:
            break
        if abs(energies_kpoint[i]-energies_kpoint[0.2]) < 0.005:
            print("This structure is converged at a k-point spacing of",i,"1/AA.")
            conv_k = xmls[i]
            break
    # exception if no convergence was reached 
    if abs(energies_kpoint[0.3]-energies_kpoint[0.2]) > 0.005:
        conv_k = xmls[0.2]
        print("KSPACING for this structure could not be converged with the given set of values. Proceeding with a spacing of 0.2 1/AA. Please check and do some additional testing.")

    # re-optimizes initial structure based on new convergence criteria, writes proper KPOINTS file
    input_copy(name,workdir+"/CONVOPT")
    os.chdir(workdir+"/CONVOPT")
    incar = vasp_i.Incar.from_file("INCAR")
    conv_k.kpoints.write_file("KPOINTS")
    incar.write_file("INCAR")
        
    print("Performing optimization with converged parameters.")
    vsub(name,workdir+"/CONVOPT")
    os.chdir(workdir)
    print("Finished converging KSPACING! Resulting optimization in CONVOPT. The k-points are:")
    print(conv_k.kpoints)
