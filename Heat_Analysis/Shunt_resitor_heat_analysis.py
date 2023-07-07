# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 23:25:34 2023

@author: rdiazpacheco
Shunt Resistor
"""
#%% Dependencies
import os
os.chdir("G:\My Drive\Code\General_Code\Heat_Analysis")
from Heat_analysis_dictionary import *
from Heat_analysis_functions import *

#%%
import numpy as np

max_current = 5000
current_ramp = 50 #A/s
time_scale = np.linspace(0,(max_current/current_ramp),200) #time
current_steps = time_scale*current_ramp

mass = 1 #kg
resistance = 0.001 #ohms
electric_power = (current_steps**2)*resistance

electric_heat = []
for i in range(0,len(time_scale)):
    electric_heat.append(electric_power[i]*((max_current/current_ramp)/200)+electric_heat[i-1])