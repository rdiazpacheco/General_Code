# -*- coding: utf-8 -*-
"""
Queen Snake Analysis
Created on Wed Nov 23 17:25:30 2022

@author: rdiazpacheco
"""

#%% Dependencies


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
from scipy import signal
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import display

#% Functions

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def Extract_voltages_one_file(filename,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue):
    li = []
    df = pd.read_csv(filename, header=header_no, skiprows = range(7,24))
    fname1 = (all_files[j].partition('\\')[2])
    fname = fname1[:-4]
    I_val = df.loc[:,I_col_name];
    li.append(I_val)
    li[0].rename("I_"+fname,inplace=True)
    for i in range(1,(1+no_ch)):
        VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
        li.append(VTap)
        li[i].rename("V"+str(i),inplace=True)            
    frame = pd.concat(li, axis=1, ignore_index= False)
    return frame

def offset_voltage_perCh(data):
    I_start = np.where(data.iloc[:,0]>0)
    I_start_ind = min(I_start[0])
    Offset_values = []
    for i in range(0,tot_ch):
        Offset_values.append(np.mean(data.iloc[0:10*(I_start_ind-10),(i+1)]))
    return Offset_values

def find_start_end_ramp_onefile(data,Istart):
    max_Iall = []
    I_ind = []
    I_str_stp = np.zeros(2)
    I_nparray = np.where(data.iloc[:,0]>Istart)
    I_ind.append(I_nparray[0])
    I_str_stp[0] = int(min(I_ind[0]))
    I_max = max(data.iloc[:,0])
    max_Ip = np.where(data.iloc[:,0] == max(data.iloc[:,0]))
    max_Iall.append(max_Ip[0])
    I_str_stp[1] = int(min(max_Iall[0]))
    return I_str_stp, I_max

def I_idx(data,I_value):
    I_ind = []
    I_nparray = np.where(data.iloc[:,0]>I_value)
    I_ind.append(I_nparray[0])
    I_indx = int(min(I_ind[0]))
    return I_indx

def range_between_two_Ivalues(data,Istart,Iend):
    I_nA = []
    I_nB = []   
    I_ABs = np.zeros(2)
    I_nAp = np.where(data.iloc[:,0] > Istart)
    I_nA.append(I_nAp[0])   
    I_ABs[0] = int(min(I_nA[0]))
    #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
    if Iend >= max(data.iloc[:,0]):  
        I_nB.append(min(I_nA[0]) + 500)
        I_ABs[1] = int(min(I_nA[0]) + 500)
    else:
        I_nBp = np.where(data.iloc[:,0] > Iend)
        I_nB.append(I_nBp[0])    
        I_ABs[1] = int(min(I_nB[0]))
    return I_ABs

def average_in_range(indices,data,filt_before):  
    Avg_noise1 = []
    if filt_before == 1:    
        for i in range(0,tot_ch):
            Avg_noise1.append(np.mean(
                butter_lowpass_filter(cutoff, fs, order,data.iloc[int(10*indices[0]):int(10*indices[1]),i+1])))           
    else:
        for i in range(0,tot_ch):
            Avg_noise1.append(np.mean(data.iloc[int(10*indices[0]):int(10*indices[1]),i+1]))
    Avg_noise = pd.Series(Avg_noise1)                  
    return Avg_noise

def average_value_at_step(data,step_currents,decay_time):
    #step_currents is an array with the step values
    aaa = []
    for i in range(0,len(step_currents)):
        range_indices = range_between_two_Ivalues(data,step_currents[i]-1,step_currents[i]+1)
        average_in_step = []
        
        for j in range(1,tot_ch+1):
            average_in_step.append(np.mean(data.iloc[int(10*(range_indices[0]+decay_time)):int(10*range_indices[1]),j]))
        aaa.append(pd.Series(average_in_step))
        average_Vstep_per_tap = pd.concat(aaa, axis=1, ignore_index= False)  
    return average_Vstep_per_tap

def func(x, Vc, Ic, V_floor, n):
    return Vc*(x/Ic)**n + V_floor

def func_w_R(x, Vc, Ic, V_floor, n,R):
    return Vc*(x/Ic)**n + V_floor + Vc*x*R


def func_only_R(x,R,OffSet):
    return x*R+OffSet

def assign_names(n1, n2, n3, n4):
    vtap_names = []
    vtap_names.append(str(n1))
    vtap_names.append(str(n2))
    vtap_names.append(str(n3))
    vtap_names.append(str(n4))
    return vtap_names
#%% Data import
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
num_files = len(all_files)
All_files = {}
no_ch = 8
for j in range(0,num_files):  
    one_filename = all_files[j]
    #All_files[str(one_filename[-25:-4])+"{0}".format(j)] = Extract_voltages_one_file(one_filename,8,4,"CmdAmps.Value","ChStat","Value")
    All_files[str(all_files[j].partition('\\')[2][-5:-4])] = Extract_voltages_one_file(one_filename,no_ch,4,"CmdAmps.Value","ChStat","Value")
    #All_files[str(j+1)].iloc[:,0].interpolate("linear", inplace = True)


#butter filter specs
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

tot_ch = no_ch


#%% Fitting without steps all in the file

Tap_dist = {
    1:21,
    2:20,
    3:6,

    4:9,
    5:12,
    6:21.5,
    7:20,
    8:21,    
    }

Tap_name = {
    1:"Lead QS_B2 Good",
    2:"QS_B2 Clamp",
    3:"QS_B2 Close",
    4:"QS_B2 Mid",
    5:"QS_B2 Far",
    6:"QS_B1-B2 PV Joint",
    7:"QS_B1",
    8:"Lead QS_B1",    
    }

F_start = int(all_files[0].partition('\\')[2][-5])
ch_no = 1
Inv_taps = 1
#Mag_f = (Mag_factor_correction[ch_no]/(1e+6))
Mag_f = 1
#Starting current
I_start = 200
#End of resisitve Current Value
num_files = len(all_files)
#ch_no = 6
R_ind = 2200
#End of noise Current Value
N_ind = 300

first_guess = [1500, 1e-6, 10, 100e-9]

for j in range(0, len(all_files)):
    File_num = j+F_start
    #Select noice range, give it a intial and final current where the noise should be
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
    #Select a range where the voltage can be considered as noise
    I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
    
    Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
    Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)
    Avg_inductive_V = Inv_taps*Mag_f*(Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1])    
    x_data = All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = signal.decimate(Inv_taps*Mag_f*(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),10)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)

    first_guess = [2000, 1e-6, 10, 100e-9]
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)

    #Display     
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.000025))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)
    
    ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0], 
            signal.decimate(Inv_taps*Mag_f*All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no],10),
            label = "Raw", linewidth = 5, linestyle = "-", color = "purple", alpha = 0.5-(j/10))
    #Raw - inductive voltage
    ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0], 
            signal.decimate(Inv_taps*Mag_f*(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),10),
            label = "Raw - V(ind)", linewidth = 2, linestyle = "-", color = "tab:green")
    
   
    #Fitted Function
    ax.plot(fit_x, fit_y,label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 0.9-(j/10))

    
    #Display crical voltage and current
    ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            np.full_like(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
            linewidth = 3, linestyle = ":", color = "red", alpha = 1, label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 18)

#%% Fitting for joint resistance: Ch 1


Tap_dist = {
    1:21,
    2:20,
    3:20,
    4:32,
    5:22.5,
    6:20,
    7:20,
    8:21,    
    }

Tap_name = {
    1:"Lead QS_B2 Good Post Epoxy",
    2:"QS_B2 Clamp",
    3:"QS_B2 Close",
    4:"QS_B2 Mid",
    5:"QS_B2 Far",
    6:"QS_B1-B2 PV Joint Post Epoxy",
    7:"QS_B1",
    8:"Lead QS_B1 Post Epoxy",    
    }
all_files.sort()
F_start = int(all_files[0].partition('\\')[2][-5])
ch_no = 8
Inv_taps = 1
Mag_f = 1
#Starting current
I_start = 300
#End of resisitve Current Value
num_files = len(all_files)

#ch_no = 6
R_ind = 1200
#End of noise Current Value
N_ind = 500

first_guess = [40e-9, 1e-5]

for j in range(0, len(all_files)):
    File_num = j+F_start
    #Select noice range, give it a intial and final current where the noise should be
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
    #Select a range where the voltage can be considered as noise
    I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
    Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)

    Avg_inductive_V = Inv_taps*Mag_f*(Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1])    
    x_data = All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = signal.decimate(Inv_taps*Mag_f*(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),10)    
   
    popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_y = func_only_R(fit_x, *popt)
    
    fig = plt.figure(1)
    
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])
    fname = fname1[:-4]
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
       
    #"""
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)
    #"""
    #Plots 
    ax.plot(x_data,y_data, linewidth = 2, linestyle = "-", color = "black", alpha = 0.5-(j/10))
    ax.plot(fit_x,fit_y,
            label = Tap_name[ch_no] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
            linewidth = 5, linestyle = "--", color = "tab:blue")
    
    ax.legend(fontsize = 25)


#%%
# Ramp information
Steps = {
    1:[3500,3700,3900,4000,4750],
    2:[3500,3700,3900,4000,4750],
    3:[4000,4200,4300,4500,4750],
    4:[4000,4200,4300,4500,4750],
    5:[4000,4200,4300,4500,4750],
    6:[4000,4150,4300,4450,4750],
    7:[4000,4150,4300,4450,5000],    
    8:[4000,4150,4300,4450,5000]
    }

Steps2 = {
    #8:[4100,4200,4300,4500,4750],
    9:[3500,3750,4000,4250,4500],
    0:[4100,4200,4300,4400,4500],
    2:[4100,4200,4300,4400,4750],
    3:[4100,4200,4300,4400,4750],
    }

Tap_dist = {
    1:120,
    2:120,
    3:20,
    4:1180,
    5:790,
    6:790,
    7:790,
    8:20,    
    }

#%% All files extracted, in a large dictionary, where each file has its own matrix


#%%
# Ramp information
Steps = {
    1:[3500,3700,3900,4000,4750],
    2:[3500,3700,3900,4000,4750],
    3:[4000,4200,4300,4500,4750],
    4:[4000,4200,4300,4500,4750],
    5:[4000,4200,4300,4500,4750],
    6:[4000,4150,4300,4450,4750],
    7:[4000,4150,4300,4450,5000],    
    8:[4000,4150,4300,4450,5000]
    }

Steps2 = {
    #8:[4100,4200,4300,4500,4750],
    9:[3500,3750,4000,4250,4500],
    0:[4100,4200,4300,4400,4500],
    2:[4100,4200,4300,4400,4750],
    3:[4100,4200,4300,4400,4750],
    }

Tap_dist = {
    1:120,
    2:120,
    3:20,
    4:1180,
    5:790,
    6:790,
    7:790,
    8:20,    
    }

tot_ch = no_ch

#%% One Channel at time stepped function
#File number
File_num = 8
#Starting current
I_start = 500
#End of resisitve Current Value
num_files = len(all_files)

#ch_no = 6
R_ind = 2000
#End of noise Current Value
N_ind = 2000

#Select noice range, give it a intial and final current where the noise should be
Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
#Select a range where the voltage can be considered as noise
I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
#Find the average noise/signal between the given ranges

ch_no = 1
Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)
# Stepped portion 
Average_signal_at_steps = average_value_at_step(All_files[str(File_num)],Steps[File_num],100)
Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1]

first_guess = [(Tap_dist[ch_no]*(1e-6)), 4300, Avg_inductive_V, 20] # First guess

#x-data here
x_data = Steps[File_num].copy()
I_rep = 4750
x_data.insert(0,500)
x_data.insert(1,N_ind)
#===============> if some points noisy
del x_data[4:6] #===============> if some points noisy

#==============> If the channel is railed
#I_rep = 4750 #==============> If the channel is railed
#x_data = np.array(x_data) #==============> If the channel is railed
#x_data[-1] = I_rep #==============> If the channel is railed

#Y_data here 
y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
y_data = list(map(float, y_data))
y_data.insert(0,All_files[str(File_num)].iloc[int(10*I_idx(All_files[str(File_num)],500)),ch_no]-Avg_inductive_V)
y_data.insert(1,All_files[str(File_num)].iloc[int(10*I_idx(All_files[str(File_num)],N_ind)),ch_no]-Avg_inductive_V)

del y_data[4:6] #===============> if some points noisy

#if the chanel is railed:
#y_data = np.array(y_data)
#I_replacement = All_files[str(File_num)].iloc[10*I_idx(All_files[str(File_num)],I_rep),ch_no] - Avg_inductive_V
#y_data[-1] = I_replacement


#fit function
popt, pcov = curve_fit(func,x_data,y_data,p0=first_guess)
fit_x = np.linspace(I_start,Imax,len(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]))
fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]))
fit_y = func(fit_x, *popt)

plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()
all_files.sort()
fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

fname = fname1[:-4]
plt.rcParams["figure.figsize"] = [25, 15]
fig.suptitle("Ch. " + str(ch_no) + " - " + fname,fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_axisbelow(True)
"""
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
#"""
ax.set_xlabel("Current [A]", fontsize=25)
ax.set_ylabel("Voltage [V]", fontsize=25)

ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0], 
        signal.decimate(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no],10),
        label = "Raw", linewidth = 5, linestyle = "-", color = "black")
#Raw - inductive voltage
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0], 
        signal.decimate(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1],10),
        label = "Raw - V(ind)", linewidth = 2, linestyle = "-", color = "tab:blue")
#Fitted Function
ax.plot(fit_x, fit_y,label = "Fit, n=" + str(round(popt[3],1)) + " Ic=" + str(round(popt[1],1)), linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 0.7)

#Display crical voltage and current
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],
        np.full_like(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],popt[0]),
        linewidth = 3, linestyle = "-.", color = "red", alpha = 0.7, label = "Fit Vc = " + str(round(popt[0]*1e6,1)) + "uV")
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],
        np.full_like(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1, label = "Real Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
ax.axvline(x = round(popt[1],1), color = 'red', linewidth = 2,linestyle='-.')
#Average at stepped locations
ax.scatter(x_data,y_data, s = 700, marker = "v", color = "orange", label = "Fit data", zorder = 4)
ax.legend(fontsize = 20)


#%% Fitting for joint resistance: Ch 1
#%% Fitting R stepped
File_num = 3
ch_no = 1
I_start = 200
R_ind = 2000
I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
#Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
#Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)
#Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1]

def func_only_R(x,R,OffSet):
    return x*R+OffSet
first_guess = [40e-9, 1e-5]



x_data = (All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0])
y_data = signal.decimate(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no],10)

popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)

fit_x = np.linspace(I_start,R_ind,len(x_data))
fit_y = func_only_R(fit_x, *popt)

all_files.sort()
fname1 = (all_files[File_num-1].partition('\\')[2])
fname = fname1[:-4]
#plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()
fig.suptitle("Ch. " + str(ch_no) + " - " + fname,fontsize=25)

"""
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.yaxis.offsetText.set_fontsize(20)
ax.set_axisbelow(True)
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(0.00001))
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.set_xlabel("Current [A]", fontsize=25)
ax.set_ylabel("Voltage [uV]", fontsize=25)
"""
#Plots 
ax.plot(x_data,y_data, linewidth = 2, linestyle = "-", color = "black")
ax.plot(fit_x,fit_y,
        label = "R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
        linewidth = 5, linestyle = "--", color = "tab:red")
ax.legend(fontsize = 40)

#%
fname = fname + " Ch"+str(ch_no)
plt.savefig(fname + ".png", dpi = 300)



