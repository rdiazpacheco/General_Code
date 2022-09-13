# -*- coding: utf-8 -*-
"""
Jabba Parse V2
Created on Wed Sep  7 15:43:57 2022

@author: rdiazpacheco

Output:
    One graph that shows:
        1. Scatter of all input in V or V/cm???
        2. A LPF of the average 
        3. A fit to the average
        4. Display N value and Ic

Give:
    1. A folder with the Ic ramps that reach Ic

This Jabba Parser will do the following:
    1. Extract data from a folder
    2. Find the beginning and end of the ramp
    3. Select a Tap, copy all data to new df
    4. Average all of them 
    
"""
#%% Functions 

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
import scipy.interpolate as spi
from matplotlib.widgets import TextBox
from matplotlib.pyplot import cm
from scipy.signal import lfilter
from scipy.signal import butter,filtfilt
from scipy.stats import linregress
import scipy as scipy
from scipy.optimize import curve_fit
import math
import sys, os, argparse
from lmfit import minimize, Parameters, Parameter, report_fit, Model
import peakutils
from scipy.interpolate import interp1d
#import xlsxwriter
import pdb # debugging tool
import re


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

def Extract_Voltages2(all_files):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):        
        df = pd.read_csv(all_files[j], header=6, skiprows = range(7,24))
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

def get_average_of_tap(tap_of_interest):
    li = [] 
    li.append(pd.Series(InVs_1.iloc[:,0].iloc[int(Current_indices[0,0]):int(Current_indices[0,1])].values))

    li1 = []
    for j in range(0,len(all_files)):
        li1.append(InVs_1.iloc[:,5*j+(tap_of_interest)].iloc[int(Current_indices[j,0]):int(Current_indices[j,1])])
        li2 = np.mean(li1, axis = 0)
    li.append(pd.Series(li2))
       
    avg_V_tap = pd.concat(li, axis=1, ignore_index= False)    
    
    return avg_V_tap

def find_start_end_ramp(Istart):
    I_ind = []
    max_Iall = []
    I_str_stp = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nparray = np.where(InVs_1.iloc[:,5*i]>Istart)
        I_ind.append(I_nparray[0])
        I_max = max(InVs_1.iloc[:,5*i])
        max_Ip = np.where(InVs_1.iloc[:,5*i] == max(InVs_1.iloc[:,5*i]))
        max_Iall.append(max_Ip[0])
        I_str_stp[i,0] = min(I_ind[i])
        I_str_stp[i,1] = min(max_Iall[i])        
    return I_str_stp, I_max

def noise_range(Istart,Iend):
    I_nA = []
    I_nB = []   

    I_ABs = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nAp = np.where(InVs_1.iloc[:,5*i] > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[i,0] = min(I_nA[i])
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(InVs_1.iloc[:,5*i]):  
            I_nB.append(min(I_nA[i]) + 500)
            I_ABs[i,1] = min(I_nA[i]) + 500
        else:
            I_nBp = np.where(InVs_1.iloc[:,5*i] > Iend)
            I_nB.append(I_nBp[0])    
            I_ABs[i,1] = min(I_nB[i])          
    return I_ABs

def find_noise_per_tap(I_ABs,filt_before):
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
    
    Avg_noise1 = []
    if filt_before == 1:    
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise1.append(np.mean(butter_lowpass_filter(cutoff, fs, order,InVs_1.iloc[int(I_ABs[j,0]):int(I_ABs[j,1]),j*5+i+1])))   
    else:
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise1.append(np.mean(InVs_1.iloc[int(I_ABs[j,0]):int(I_ABs[j,1]),j*5+i+1]))
    Avg_noise = pd.Series(Avg_noise1)         
    
    avg_noise_p_tap = []           
    for i in range(0,len(all_files)):
        avg_noise_p_tap.append(np.mean(Avg_noise[[i,i+4,i+8,i+12]]))
      
    return Avg_noise, avg_noise_p_tap

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def func(x, Vc, Ic, V_floor, n):
    return Vc*(x/Ic)**n + V_floor

def assign_names(n1, n2, n3, n4):
    vtap_names = []
    vtap_names.append(str(n1))
    vtap_names.append(str(n2))
    vtap_names.append(str(n3))
    vtap_names.append(str(n4))
    return vtap_names

def find_Ic_values(a,b,c,d,avg_npt):
           
    Ic_values = np.zeros([4*len(all_files),2])
    V_over_1 = []
    for j in range(0,len(all_files)):
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        
        for i in range(0,4):
            #V_over1 = np.where(InVs_1["VTap"+str(i+1)+"_"+fname]>1)
            V_over1 = np.where(butter_lowpass_filter(cutoff, fs, order,
                                                     (InVs_1["VTap"+str(i+1)+"_"+fname].iloc[0:int(Current_indices[j,1])]))-avg_npt[(4*j)+i]>1)
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
    Ic_avg = []
    for i in range(0,len(all_files)):
        Ic_avg.append(np.mean(Ic_values[[i,i+4,i+8,i+12],1]))
    
    return Ic_values, Ic_avg

def func(x, Vc, Ic, V_floor, n):
    return Vc*(x/Ic)**n + V_floor
#%%
#####
####
#Here you can visualize all the taps in all the files. It will open as many pictures as there a files. Run this section to pull the data.
# Then run the section below to visualize 
cutoff = 1     # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 1     # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples
#extract all data from a folder
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
InVs_1 = Extract_Voltages(all_files)

#Select noice range
Noice_ind = noise_range(200, 1000)
Avg_noises, n_p_t_avg = find_noise_per_tap(Noice_ind,1)
Current_indices, Imax = find_start_end_ramp(200)
Ic_vals, Ic_avg = find_Ic_values(0,0,1,0,Avg_noises)



#%%
# Run this section to visualize all your taps

# give names to the taps
vtns = assign_names("Cable: A1", "Cable: A2", "PIT-V Joint", "Overall") #for 20220819
 

for k in range(0,len(all_files)):
#for k in range(0,4):
    fig, ax = plt.subplots(2,2, figsize=(25, 35))

    fname1 = (all_files[k].partition('\\')[2])
    fname = fname1[:-4]
    #fig.suptitle('Jumper Cables Heat Analysis',fontsize=25)
    fig.suptitle(fname,fontsize=25)
    
 
    
    SC_parameter = np.zeros(len(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]))
    SC_parameter[:] = 1

    for j in range(0,2):
        for i in range(0,2):
            ax[i,j].set_title("Channel "+str(i+1+(2*j))+" - "+ vtns[i+(2*j)],fontsize=20)
            ax[i,j].set_xlabel("Current [A]", fontsize=15)
            ax[i,j].set_ylabel("Electric Field [uV/cm]", fontsize=15)
            #ax[i,j].set_xscale('log')
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
                       InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*k)+(2*j)+i],
                       label="Voltage Tap", linewidth = 3,c = c)
            c = next(color) 
            c = next(color) 
            IV_1, = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                        butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                                              .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*k)+(2*j)+i]),
                        label="2Hz LPF", linewidth = 3,c = c)
            
            OneuV = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                                 SC_parameter, label="1uV/cm", color='red', linewidth=4, linestyle=':')
            
    
    for j in range(0,2):
        for i in range(0,2):           
            leg = ax[i,j].legend(fontsize=15)

##% extract all data from a folder
# This section is the same as above

sample_name = "Glycon II: "
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
InVs_1 = Extract_Voltages(all_files)
#For Glycon II
#InVs_1 = Extract_Voltages2(all_files)
#%%
# This will focus on one tap and do all the fittings and find Ic and N values (the fitting)
# 1. Select tap number (1,2,3,4)
# Select a guess for the n value, usually something between 5-10. You can reiterate to get a better fit. If it's too far (over or under, the fit will be obviosly bad)
# give the sample a name "queen snake"
# and that's it. It will average all the runs of the same kind of the same filter and do all the analysis.

#Select noice range
sample_name = "Queen Snake: "
tap_of_interest = 4
n_guess =5

#butter low pass filter is used here
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 1     # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

Noice_ind = noise_range(200, 1000)
Avg_noises, n_p_t_avg = find_noise_per_tap(Noice_ind,1)
Current_indices, Imax = find_start_end_ramp(200)
Ic_vals, Ic_avg = find_Ic_values(1,1,1,1,Avg_noises)

#get average run
Avg_InV_tap = get_average_of_tap(tap_of_interest)

#first guess 
first_guess = [1, Ic_avg[int(tap_of_interest)-1], n_p_t_avg[int(tap_of_interest)-1],n_guess]
#first_guess = [1, 2500, n_p_t_avg[int(tap_of_interest)-1],n_guess]
#fit function
#popt, pcov = curve_fit(func, Avg_InV_tap.iloc[:,0], Avg_InV_tap.iloc[:,1], p0=first_guess)
popt, pcov = curve_fit(func, butter_lowpass_filter(cutoff, fs, order,Avg_InV_tap.iloc[:,0]), Avg_InV_tap.iloc[:,1], p0=first_guess)

fit_x = np.linspace(200,Imax,3000)
fit_y = func(fit_x, *popt)

SC_parameter = np.zeros(len(Avg_InV_tap))
SC_parameter[:] = 1                     

plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()

fname1 = (all_files[i].partition('\\')[2])
fname = fname1[:-4]
fig.suptitle(sample_name + fname[:-4],fontsize=25)

ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

#plt.figure(figsize = (30,18))

ax.set_ylim((-1,80))
ax.set_xlim((50,int(Imax+50)))
ax.set_xlabel("Current [A]", fontsize=15)
ax.set_ylabel("Electric Field [uV/cm]", fontsize=15)



ax.plot(Avg_InV_tap.iloc[:,0],butter_lowpass_filter(cutoff, fs, order,Avg_InV_tap.iloc[:,1])-n_p_t_avg[int(tap_of_interest)-1], 
         label = "Average Raw", linewidth = 5, linestyle = "-", color = "tab:blue")
ax.plot(fit_x, fit_y-n_p_t_avg[int(tap_of_interest)-1],
         label = "Fit, n=" + str(round(popt[3],1)) + " Ic=" + str(round(popt[1],1)), linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 0.7)
for i in range(0,len(all_files)):
    fname1 = (all_files[i].partition('\\')[2])
    fname = fname1[:-4]
    ax.scatter(InVs_1.iloc[int(Current_indices[i,0]):int(Current_indices[i,1]),5*i],
                           InVs_1.iloc[int(Current_indices[i,0]):int(Current_indices[i,1]),5*i+tap_of_interest]-n_p_t_avg[int(tap_of_interest)-1],
                           label = "Raw " + fname,marker= 10, s= 6, alpha=0.2)


ax.plot(Avg_InV_tap.iloc[:,0],SC_parameter, label="1uV/cm", color='red', linewidth=2, linestyle=':')
ax.axvline(x = round(popt[1],1), color = 'red', linewidth = 2,linestyle=':')

ax.legend(fontsize = 20)

fig.show()

#%%

    
