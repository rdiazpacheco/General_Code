# -*- coding: utf-8 -*-
"""
Jabba Ic Analyzer for Glycon 1-4
Created on Thu Jan 12 12:30:53 2023

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
    try:
        I_val = df.loc[:,I_col_name];
    except: 
        header_no = 4
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

def offset_voltage_perCh(data,IvsV):
    I_start = np.where(data.iloc[:,0]>0)
    I_start_ind = min(I_start[0])
    Offset_values = []
    for i in range(0,tot_ch):
        Offset_values.append(np.mean(data.iloc[0:IvsV*(I_start_ind-10),(i+1)].astype(float)))
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

def average_in_range(indices,data,filt_before,IvsV):  
    Avg_noise1 = []
    if filt_before == 1:    
        for i in range(0,tot_ch):
            Avg_noise1.append(np.mean(
                butter_lowpass_filter(cutoff, fs, order,data.iloc[int(IvsV*indices[0]):int(IvsV*indices[1]),i+1].astype(float))))           
    else:
        for i in range(0,tot_ch):
            Avg_noise1.append(np.mean(data.iloc[int(IvsV*indices[0]):int(IvsV*indices[1]),i+1].astype(float)))
    Avg_noise = pd.Series(Avg_noise1)                  
    return Avg_noise

def average_value_at_step(data,step_currents,decay_time,IvsV):
    #step_currents is an array with the step values
    aaa = []
    for i in range(0,len(step_currents)):
        range_indices = range_between_two_Ivalues(data,step_currents[i]-1,step_currents[i]+1)
        average_in_step = []
        
        for j in range(1,tot_ch+1):
            average_in_step.append(np.mean(data.iloc[int(IvsV*(range_indices[0]+decay_time)):int(IvsV*range_indices[1]),j]))
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


#%% Glycon I
Tap_dist = {
    1:40,
    2:70,
    3:22,
    4:40    
    }
#"""


#""" 20221116
Tap_name = {
    1:"Glycon I: Lead A",
    2:"Glycon I: Core ",
    3:"Glycon I: Middle Core",
    4:"Glycon I: Lead B"
    }
#"""

Ic_RG =4130

#%% Glycon II
Tap_dist = {
    1:44,
    2:22,
    3:59,
    4:43    
    }
#"""


#""" 20221116
Tap_name = {
    1:"Glycon II: Lead A",
    2:"Glycon II: Compressed Core ",
    3:"Glycon II: Whole Core",
    4:"Glycon II: Lead B"
    }
#"""
Ic_RG =4130
#%% Glycon III
Tap_dist = {
    1:22.75,
    2:24.95,
    3:47.9,
    4:20    
    }
#"""


#""" 20221116
Tap_name = {
    1:"QS_A1-QS_A2 PVJ 1",
    2:"QS_A1-QS_A2 PVJ 2",
    3:"QS_A1_479",
    4:"QS_A1_200"
    }
#"""

#%% Glycon IV
Tap_dist = {
    1:22.75,
    2:24.95,
    3:47.9,
    4:20    
    }
#"""


#""" 20221116
Tap_name = {
    1:"QS_A1-QS_A2 PVJ 1",
    2:"QS_A1-QS_A2 PVJ 2",
    3:"QS_A1_479",
    4:"QS_A1_200"
    }
#"""

#%% Data Import
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
num_files = len(all_files)
All_files = {}
no_ch = 4
header_no = 6

for j in range(0,num_files):  
    one_filename = all_files[j]
    All_files[str(all_files[j].partition('\\')[2][-5:-4])] = Extract_voltages_one_file(one_filename,no_ch,header_no,"CmdAmps.Value","VTapFilter","rVal")
tot_ch = no_ch



#%%
all_files.sort()
F_start = int(all_files[0].partition('\\')[2][-5])
ch_no = 3
Inv_taps = 1
#Mag_f = (Mag_factor_correction[ch_no]/(1e+6))
Mag_f = 1e-7
#Starting current
I_start = 150
#End of resisitve Current Value
num_files = len(all_files)
#ch_no = 6
R_ind = 4599
#End of noise Current Value
N_ind = 4199
IvsV = 1

first_guess = [3500, 1e-6, 11, 50e-9]
ax.cla()
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[str(File_num)],0,1)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[str(File_num)],0,1)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[str(File_num)],1)[ch_no-1])
    
    x_data = All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*(All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)

    first_guess = [2400, 1e-6, 10, 100e-9]
    
    #Display     
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=30)
    ax.set_ylabel("Voltage [V]", fontsize=30)
    #Raw 
    #ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0], 
    #        signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no],IvsV),
    #        label = "Raw", linewidth = 0.1, linestyle = "-", color = "tab:blue", alpha = 0.75-(j/10))
    #Cleaned - inductive voltage
    ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0], 
            signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*
                            (All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV),
            label = "Raw - V(ind)", linewidth = 0.25, alpha = 0.85-(j/10), linestyle = "-", color = "black")
    
   
    #Fitted Function
    ax.plot(fit_x, fit_y,label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 1-(j/10))

    
    #Display crical voltage and current
    ax.plot(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            np.full_like(All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
            linewidth = 3, linestyle = ":", color = "red", alpha = 1, 
            label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV, Degradaton = " + str(round((100*(1-(popt[0]/Ic_RG))),2)) + "%")
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 15)
    
#%% Fitting for joint resistance: Ch 1
all_files.sort()
F_start = int(all_files[0].partition('\\')[2][-5])
ch_no = 4
Inv_taps = 1
Mag_f = 1e-7
#Starting current
I_start = 300
#End of resisitve Current Value
num_files = len(all_files)

#ch_no = 6
R_ind = 2500
#End of noise Current Value
N_ind = 500

first_guess = [40e-9, 1e-5]
ax.cla()
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[str(File_num)],0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[str(File_num)],0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[str(File_num)],IvsV)[ch_no-1])
    
    x_data = All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*
                             (All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]),IvsV)
   
    popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
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
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    ax.set_xlabel("Current [A]", fontsize=30)
    ax.set_ylabel("Voltage [V]", fontsize=30)
    #"""
    #Plots 
    ax.plot(x_data,y_data, linewidth = 0.25, linestyle = "-", color = "black", alpha = 0.95-(j/10))
    ax.plot(fit_x,fit_y,
            label = Tap_name[ch_no] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
            linewidth = 5, linestyle = "--", color = "tab:red")
    
    ax.legend(fontsize = 40)





