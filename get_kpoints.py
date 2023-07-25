#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 16:03:47 2023

@author: bt308570
"""

# simple script to extract the k-point settings from the vasprun.xml file

from pymatgen.io.vasp.outputs import Vasprun
import sys

name = sys.argv[1]
xml = Vasprun(name)
print(xml.kpoints.as_dict)