# -*- coding: utf-8 -*-
"""
Transverse Compression Test 
Critical Current  finder
    - CSV file parser
    - Uses Critical crrent code from Brandon and Zach and Erica with some modifications
"""
#%% Preliminaries

import csv
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.interpolate as spi
from matplotlib.widgets import TextBox
import os
from os import listdir
from os.path import isfile, join



#%%Data location 
#print(os.getcwd())
filespath = 'CFS Dropbox\\CFS_Internal\\R&D\\CSMC\\Testing\\CSMC-176 Transverse Compression Testing\\04 Test Data'
os.chdir(filespath)  # Provide the new path here
#%%
print(os.getcwd())

#Obtain list of available files
filelist = [f for f in listdir('\Glycon_III') if isfile(join(filespath, f))]

#%% CSV Parser

with open('employee_birthday.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            line_count += 1
    print(f'Processed {line_count} lines.')
