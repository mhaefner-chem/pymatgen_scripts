#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 11:55:16 2023

@author: bt308570
"""

import subprocess

process = subprocess.Popen(['squeue',"-u","bt308570"],
             stdin=subprocess.DEVNULL,
             stdout=subprocess.PIPE,
             start_new_session=False)
out, err = process.communicate()
jobs = out.decode().count("bt308570")

print(jobs)