# -*- coding: utf-8 -*-
"""
Jabba the Cart Parse for QMC hetaers
Based on general Jabba Parser
Created on Sat Aug 27 07:08:20 2022

@author: rdiazpacheco
"""
import csv
import time
import json
import copy
import pandas as pd
from tkinter import filedialog
from influxdb import InfluxDBClient
import glob
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.interpolate as spi
from matplotlib.widgets import TextBox
from matplotlib.pyplot import cm
from scipy.signal import lfilter
from scipy.signal import butter,filtfilt

def Extract_Voltages_only(all_files,No_of_taps):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):
        
        df = pd.read_csv(all_files[j], header=4, skiprows = range(7,24))
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        #I_val = df.loc[:,"CmdAmps.Value"];
        #li.append(I_val)
        #li[(5*j)].rename("I_"+fname,inplace=True)
        if No_of_taps == 1:
            VTap = df.loc[:,"VTapFilter["+str(1)+"].rVal"];
            li.append(VTap)
                  
        else:
            for i in range(1,No_of_taps):
                VTap = df.loc[:,"VTapFilter["+str(i)+"].rVal"];
                li.append(VTap)
                li[(i)+(No_of_taps+1*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)
                
        li[j].rename("VTap"+str(1)+"_"+fname,inplace=True) 
        frame = pd.concat(li, axis=1, ignore_index= False)        
    return frame

def find_start_end_ramp2(No_of_taps,n_jump,v_jump,InVs_1):
    
    Volt_jumps = []
    vjump = pd.DataFrame()
    for j in range(0,len(all_files)):
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        if No_of_taps == 1:
            Volt_jumps1 = np.zeros([(int(len(Vs_1["VTap"+str(1)+"_"+fname])/n_jump)),1])
            for i in range(0,int(len(Vs_1["VTap"+str(1)+"_"+fname])/n_jump)):
                Volt_jumps1[i,0] = Vs_1["VTap"+str(1)+"_"+fname].iloc[n_jump*i+4]-Vs_1["VTap"+str(1)+"_"+fname].iloc[n_jump*i]
            Volt_jumps.append(Volt_jumps1)
                           
        else: 
            Volt_jumps1 = np.zeros(len(Vs_1.loc["VTap"+str(1)+"_"+fname]),1)
            for k in range(0,No_of_taps):
                
                for i in range(0,len(Vs_1.loc["VTap"+str(1)+"_"+fname])):
                    Volt_jumps1[i,j] = Vs_1["VTap"+str(k)+"_"+fname].iloc[5*i+4]-Vs_1["VTap"+str(k)+"_"+fname].iloc[5*i]
            Volt_jumps.append(Volt_jumps1)

        
        for k in range(0,len(Volt_jumps)):
            vjump['filename'] = np.where(abs(Volt_jumps[k])>v_jump)
            
         
    return vjump

#%%
#%% Extract all data from a day
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
#%%

Vs_1 = Extract_Voltages_only(all_files,1) #extract all data
ramps = find_start_end_ramp2(1,5,5,Vs_1) #find the places where the ramp begins and ends

#%% save as coordinates
No_of_taps = 1
n_jump = 5
v_jump = 7
Volt_jumps = []
vjump2 = []

for j in range(0,len(all_files)):
    fname1 = (all_files[j].partition('\\')[2])
    fname = fname1[:-4]
    if No_of_taps == 1:
        Volt_jumps1 = np.zeros([(int(len(Vs_1["VTap"+str(1)+"_"+fname])/n_jump)),1])
        for i in range(0,int(len(Vs_1["VTap"+str(1)+"_"+fname])/n_jump)):
            Volt_jumps1[i,0] = Vs_1["VTap"+str(1)+"_"+fname].iloc[n_jump*i+4]-Vs_1["VTap"+str(1)+"_"+fname].iloc[n_jump*i]
        Volt_jumps.append(Volt_jumps1)
                       
    else: 
        Volt_jumps1 = np.zeros(len(Vs_1.loc["VTap"+str(1)+"_"+fname]),1)
        for k in range(0,No_of_taps):
            
            for i in range(0,len(Vs_1.loc["VTap"+str(1)+"_"+fname])):
                Volt_jumps1[i,j] = Vs_1["VTap"+str(k)+"_"+fname].iloc[5*i+4]-Vs_1["VTap"+str(k)+"_"+fname].iloc[5*i]
        Volt_jumps.append(Volt_jumps1)

#%%
vjump2_sizes = []
for j in range(0,len(Volt_jumps)):
    vjump1 = np.where(abs(Volt_jumps[j])>v_jump)
    vjump2.append(n_jump*vjump1[0])
    #vjump2_sizes.append(len(vjump2[j]))
    #vjump2[j].rename("VTap"+str(i)+"_"+fname,inplace=True)
    
#%%
# specifying the plot size
plt.figure(figsize = (10, 5))
 
# only one line may be specified; full height
time = np.linspace(1,len(Vs_1.iloc[:,1]),len(Vs_1.iloc[:,1]))
plt.plot(time,Vs_1.iloc[:,1],linewidth = 3)

for i in range(0,len(vjump2[1])):
    plt.axvline(x = vjump2[1][i], color = 'r', label = 'Start/end pulse',linewidth = 1, linestyle = ":")
 
# rendering plot
plt.show()

#%%
plt.figure(figsize = (20, 20))
 
# only one line may be specified; full height

offset1 = 1
offset2 = 7

for i in range(0,int(((len(vjump2[1]))/2)-2)):
    time = np.linspace(1,len(Vs_1.iloc[vjump2[1][2*i]+offset1:vjump2[1][2*i+1]+offset2,1]),((vjump2[1][2*i+1]+offset2)-(vjump2[1][2*i]+offset)))
    plt.plot(time,Vs_1.iloc[vjump2[1][2*i]+offset1:vjump2[1][2*i+1]+offset2,1],linewidth = 3)
    #plt.axvline(x = vjump2[1][i], color = 'r', label = 'Start/end pulse',linewidth = 1, linestyle = ":")
    #plt.plot(time,Vs_1.iloc[:,1],linewidth = 3)
 
# rendering plot
plt.show()

plt.figure(figsize = (10, 5))
 
# only one line may be specified; full height

offset1 = 15
offset2 = 0

for i in range(0,int(((len(vjump2[1]))/2)-2)):
    time = np.linspace(1,len(Vs_1.iloc[vjump2[1][2*i]+offset1:vjump2[1][2*i+1]+offset2,1]),((vjump2[1][2*i+1]+offset2)-(vjump2[1][2*i]+offset1)))
    plt.plot(time,Vs_1.iloc[vjump2[1][2*i]+offset1:vjump2[1][2*i+1]+offset2,1],linewidth = 3)
    #plt.axvline(x = vjump2[1][i], color = 'r', label = 'Start/end pulse',linewidth = 1, linestyle = ":")
    #plt.plot(time,Vs_1.iloc[:,1],linewidth = 3)
 
# rendering plot
plt.show()

#%% 
li = []
for i in range(0,int(((len(vjump2[1]))/2)-1)):
    li.append(Vs_1.iloc[vjump2[1][2*i]:vjump2[1][2*i+1],1])

#%%
#vjump = pd.DataFrame(index = range(max(vjump2_sizes)),columns = range(len(vjump2)))
vjump = pd.concat(li, axis=1, ignore_index= False)    


#vjump = pd.DataFrame(vjump2)

#%% Find average voltage between ramps

for i in range(0,len(ramps)):
    

#%% Graph




