# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 08:32:25 2022
Jabba CSV Output Parser

@author: rdiazpacheco
From Dan Nash
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


# get start time from string, returns something like time_ns() in UTC
def parse_t0(csv_row):
    # tstr: YEAR MONTH DAY 24HR MINUTE SECOND"
    yr = int(csv_row[4])
    mo = csv_row[3].split(' ')[-2]
    day = int(csv_row[3].split(' ')[-1])
    hr = int(csv_row[5].split(":")[0])
    min = int(csv_row[5].split(":")[1])
    sec = float(csv_row[5].split(":")[2])
    tstr = "%i %s %i %i %i %i" % (yr, mo.upper(), day, hr, min, int(sec))
    t0_loc = time.strptime(tstr, '%Y %B %d %H %M %S')
    t0_UTC = time.mktime(t0_loc)
    return 1e9*(t0_UTC + (sec % 1)) # convert to ns

# parse tag header row
def parse_tags(csv_row):
    tags = []
    for i in range(int(len(csv_row)/2)):
        tags.append(csv_row[i*2+1])
    return tags

def Extract_Voltages(all_files):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):
        
        df = pd.read_csv(all_files[j], header=4, skiprows = range(7,24))
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        I_val = df.loc[:,"CmdAmps.Value"];
        li.append(I_val)
        li[(5*j)].rename("I_"+fname,inplace=True)
        for i in range(1,5):
            VTap = df.loc[:,"VTapFilter["+str(i)+"].rVal"];
            li.append(VTap)
            li[(i)+(5*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)            
        frame = pd.concat(li, axis=1, ignore_index= False)        
    return frame

def find_start_end_ramp(Istart):
    I_ind = []
    max_Iall = []
    I_str_stp = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nparray = np.where(InVs_1.iloc[:,5*i]>Istart)
        I_ind.append(I_nparray[0])
        max_Ip = np.where(InVs_1.iloc[:,5*i] == max(InVs_1.iloc[:,5*i]))
        max_Iall.append(max_Ip[0])
        I_str_stp[i,0] = min(I_ind[i])
        I_str_stp[i,1] = min(max_Iall[i])        
    return I_str_stp

def noise_range(Istart,Iend):
    I_nA = []
    I_nB = []    
    I_ABs = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nAp = np.where(InVs_1.iloc[:,5*i] > Istart)
        I_nA.append(I_nAp[0])
        I_nBp = np.where(InVs_1.iloc[:,5*i] > Iend)
        I_nB.append(I_nBp[0])
        I_ABs[i,0] = min(I_nA[i])
        I_ABs[i,1] = min(I_nB[i])
         
    return I_ABs

def find_noise_per_tap(I_ABs):
    Avg_noise = []
    for j in range(0,len(all_files)):
        for i in range(0,4):
            Avg_noise.append(np.mean(InVs_1.iloc[int(I_ABs[j,0]):int(I_ABs[j,1]),j*5+i+1]))        
    return Avg_noise

def butter_lowpass_filter(cutoff, fs, order,data):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def find_Ic_values(a,b,c,d):
    #butter low pass filter is used here
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples

    Ic_values = np.zeros([4*len(all_files),2])
    V_over_1 = []
    for j in range(0,len(all_files)):
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        
        for i in range(0,4):
            #V_over1 = np.where(InVs_1["VTap"+str(i+1)+"_"+fname]>1)
            V_over1 = np.where(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1)+"_"+fname].iloc[0:int(Current_indices[j,1])])>1)
            if len(V_over1[0]) > 1:
                V_over_1.append(V_over1[0])
            else:
                V_over_1.append([0])
                
    for i in range(0,len(V_over_1)):
        Ic_values[i,0] = min(V_over_1[i])
       
    selected_graphs = (a,b,c,d)
    for j in range(0,len(all_files)):
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        for i in range(0,len(selected_graphs)):
            if selected_graphs[i] == 1:
                Ic_values[i+4*j,1] = InVs_1["I_"+fname].iloc[int(Ic_values[i+4*j,0])]
            else:
                Ic_values[i+4*j,1] = 0
                
    return Ic_values

#%% Extract all data from a day
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
InVs_1 = Extract_Voltages(all_files)
Current_indices = find_start_end_ramp(200)

# Filtering noise, filtering from 300 - 600 amps
Noice_indices = noise_range(300, 600)
Avg_noices = find_noise_per_tap(Noice_indices)
Ic_values = find_Ic_values(0, 0, 0, 1)

            
#%%
for k in range(0,len(all_files)):
#for k in range(0,4):
    fig, ax = plt.subplots(2,2, figsize=(25, 17))

    fname1 = (all_files[k].partition('\\')[2])
    fname = fname1[:-4]
    #fig.suptitle('Jumper Cables Heat Analysis',fontsize=25)
    fig.suptitle(fname,fontsize=25)
    
 
    
    SC_parameter = np.zeros(len(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]))
    SC_parameter[:] = 1

    for j in range(0,2):
        for i in range(0,2):
            ax[i,j].set_title("Channel "+str(i+1+2*j),fontsize=20)
            ax[i,j].set_xlabel("Current [A]", fontsize=15)
            ax[i,j].set_ylabel("Electric Field [V/cm]", fontsize=15)
            leg = ax[i,j].legend(fontsize=20)
            leg = ax[i,j].legend(fontsize=20)
    
    plt.subplots_adjust(bottom=0.1,right=0.9)
    
    color = iter(cm.rainbow(np.linspace(0, 1.5, 20))) 
    
    #butter low pass filter
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
   
    
    for j in range(0,2):
        for i in range(0,2):
            c = next(color)           
            IV_1, = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                       InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                       label="Voltage Tap", linewidth = 3,c = c)
            c = next(color) 
            c = next(color) 
            IV_1, = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                        butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]),
                        label="2Hz LPF, IC = "+ str(Ic_values[4*k+2*j+i,1]), linewidth = 3,c = c)
            
            OneuV = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                                 SC_parameter, label="1uV/cm", color='red', linewidth=4, linestyle=':')
    
    for j in range(0,2):
        for i in range(0,2):           
            leg = ax[i,j].legend(fontsize=15)
  




