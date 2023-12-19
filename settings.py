#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:42:03 2023

@author: bt308570
"""

#!/usr/bin/env python3
def incar_settings(arg):
# various INCAR settings conveniently compiled in one place
    prefunc = "PE"
    func = "PE"
    disp = "D3BJ"
    
    settings = {}
    
    # preliminary settings
    settings["ISYM"] = 1
    settings["NSW"] = 0
    settings["SIGMA"] = 0.02
    settings["ISMEAR"] = 0
    settings["POTIM"] = 0.2
    
    settings["ENCUT"] = 650
    settings["ENAUG"] = 1300
    
    
    if func == "R2SCAN" or func == "SCAN":
        settings["METAGGA"] = func
    else:
        settings["GGA"] = func
        settings["ENCUT"] = 520
        settings["ENAUG"] = 1040
    
    if func == "R2SCAN" and disp == "D4":
        settings["IVDW"] = 13
        settings["VDW_S6"] = 1.000
        settings["VDW_S8"] = 0.6018749
        settings["VDW_A1"] = 0.51559235
        settings["VDW_A3"] = 5.77342911
    elif func == "R2SCAN" and disp == "rVV10":
        settings["IVDW"] = 0
        settings["LUSE_VDW"] = ".TRUE."
        settings["BPARAM"] = 11.95
        settings["CPARAM"] = 0.0093
        settings["LASPH"] = ".TRUE."
    elif disp == "D4":
        settings["IVDW"] = 13
    elif disp == "D3BJ":
        settings["IVDW"] = 12
    elif disp == "":
        settings["IVDW"] = 0
    
    if arg == "phonopy_sps":

        
        settings["KSPACING"] = 0.3
        settings["EDIFF"] = 1e-07
        
        
    elif arg == "SCANopt":
        settings["KSPACING"] = 0.3
        settings["EDIFFG"] = -0.001
        settings["EDIFF"] = 1e-07
        
    elif arg == "PSopt":
        settings["GGA"] = prefunc
        settings["EDIFF"] = 1e-7
        settings["EDIFFG"] = -0.003
        settings["NELMIN"] = 4
        settings["SYMPREC"] = 1e-4
        
    elif arg == "Fastopt":
        settings["GGA"] = prefunc
        settings["EDIFF"] = 0.0005
        settings["EDIFFG"] = 0.05
        settings["NELMIN"] = 4
        settings.pop("ENCUT")
        settings.pop("ENAUG")
        settings["PREC"] = "Normal"
        settings["SYMPREC"] = 1e-4
    
    return settings
