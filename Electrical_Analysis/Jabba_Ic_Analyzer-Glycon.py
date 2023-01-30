# -*- coding: utf-8 -*-
"""
Jabba Ic Analyzer for Glycon 1-4
Created on Thu Jan 12 12:30:53 2023

@author: rdiazpacheco
"""
#%% Dependencies
import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from QS_Data_Index import *

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

def Extract_voltages_one_file(filename,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue,load_cell):   
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
    try:
        load_cell_val = df.loc[:,load_cell];
        load_cell_val = load_cell_val.replace(r'^\s*$', np.nan, regex=True)
        li.append(pd.Series(max(load_cell_val.astype(float))))
        li[i+1].rename("Max Load",inplace=True)
    except:
        li.append(pd.Series(0))
        li[i+1].rename("Max Load",inplace=True)        
    frame = pd.concat(li, axis=1, ignore_index= False)   
    test_dict = {}
    test_dict[0] = frame
    return frame

def offset_voltage_perCh(data,IvsV):
    data = data.replace(r'^\s*$', np.nan, regex=True)
    I_start = np.where(data.iloc[:,0].astype(float)>0)
    I_start_ind = min(I_start[0])
    Offset_values = []
    for i in range(0,tot_ch):
        Offset_values.append(np.mean(data.iloc[0:IvsV*(I_start_ind-10),(i+1)].astype(float)))
    return Offset_values



def find_start_end_ramp_onefile(data,Istart):
    max_Iall = []
    I_ind = []
    I_str_stp = np.zeros(2)
    try: 
        I_nparray = np.where(data.iloc[:,0].astype(float)>Istart)
        I_ind.append(I_nparray[0])
        I_str_stp[0] = int(min(I_ind[0]))
        I_max = max(data.iloc[:,0].astype(float))
        max_Ip = np.where(data.iloc[:,0].astype(float) == max(data.iloc[:,0].astype(float)))
        max_Iall.append(max_Ip[0])
        I_str_stp[1] = int(min(max_Iall[0]))
    except:
        data = data.replace(r'^\s*$', np.nan, regex=True)
        I_nparray = np.where(data.iloc[:,0].astype(float)>Istart)
        I_ind.append(I_nparray[0])
        I_str_stp[0] = int(min(I_ind[0]))
        I_max = max(data.iloc[:,0].astype(float))
        max_Ip = np.where(data.iloc[:,0].astype(float) == max(data.iloc[:,0].astype(float)))
        max_Iall.append(max_Ip[0])
        I_str_stp[1] = int(min(max_Iall[0]))
    return I_str_stp, I_max

def I_idx(data,I_value):
    I_ind = []
    I_nparray = np.where(data.iloc[:,0].astype(float)>I_value)
    I_ind.append(I_nparray[0])
    I_indx = int(min(I_ind[0]))
    return I_indx

def range_between_two_Ivalues(data,Istart,Iend):
    I_nA = []
    I_nB = []   
    I_ABs = np.zeros(2)
    try:
        I_nAp = np.where(data.iloc[:,0].astype(float) > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[0] = int(min(I_nA[0]))
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(data.iloc[:,0].astype(float)):  
            I_nB.append(min(I_nA[0]) + 500)
            I_ABs[1] = int(min(I_nA[0]) + 500)
        else:
            I_nBp = np.where(data.iloc[:,0].astype(float) > Iend)
            I_nB.append(I_nBp[0])    
            I_ABs[1] = int(min(I_nB[0]))
    except:
        data = data.replace(r'^\s*$', np.nan, regex=True)
        I_nAp = np.where(data.iloc[:,0].astype(float) > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[0] = int(min(I_nA[0]))
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(data.iloc[:,0].astype(float)):  
            I_nB.append(min(I_nA[0]) + 500)
            I_ABs[1] = int(min(I_nA[0]) + 500)
        else:
            I_nBp = np.where(data.iloc[:,0].astype(float) > Iend)
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
    2:59,
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
IvsV = 10
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

Cycle_count = {
    1:0,
    2:0,
    3:0,
    4:0,
    5:0,
    6:50,
    7:100,
    8:200
    }

Ic_RG =4130
IvsV = 1
Mag_f = 1e-7
#%% Glycon III
Tap_dist = {
    1:45,
    2:59,
    3:22,
    4:45    
    }
#"""


#""" 20221116
Tap_name = {
    1:"Glycon III: Lead A",
    2:"Glycon III: Whole Core ",
    3:"Glycon III: Compressed Core",
    4:"Glycon III: Lead B"
    }
#"""

Cycle_count = {
    1:0,
    2:0,
    3:0,
    4:1,
    5:50,
    6:50,
    7:100,
    8:100,
    9:200,
    10:200,
    11:200
    }

Cycle_Pressure = { 
    1:"0",
    2:"0",
    3:"0",
    4:"250",
    5:"250",
    6:"275",
    7:"275",
    8:"300",
    9:"300",
    10:"325",
    11:"325",
    12:"350",
    13:"350",
    14:"375",
    15:"375",
    16:"400",
    17:"400"    
    }

Ic_RG =4300
IvsV = 1
Mag_f = 1e-7
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
all_files.sort()
num_files = len(all_files)
All_files = {}
no_ch = 4
header_no = 6

for j in range(0,num_files):  
    one_filename = all_files[j]
    #All_files[j] = Extract_voltages_one_file(one_filename,no_ch,header_no,"FdbkAmps.Value","VTapFilter","rVal","Load Cell [MPa]")
    All_files[j] = Extract_voltages_one_file(one_filename,no_ch,header_no,"CmdAmps.Value","VTapFilter","rVal","Load Cell [MPa]")
    #All_files[str(all_files[j].partition('\\')[2][-5:-4])] = Extract_voltages_one_file(one_filename,no_ch,header_no,"CmdAmps.Value","VTapFilter","rVal","Load Cell [MPa]")
tot_ch = no_ch

#FdbkAmps.Value

#%%
all_files.sort()
#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 2
Inv_taps = 1
#Mag_f = (Mag_factor_correction[ch_no]/(1e+6))
#Starting current
I_start = 200
num_files = len(all_files)
#ch_no = 6
R_ind = 4600
#End of noise Current Value
N_ind = 2000
Ic_P = np.zeros([len(all_files),2])

first_guess = [3500, 1e-6, 11, 50e-9]
try: 
    ax.cla()
except:
    j = 0
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[int(File_num)],0,1)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[int(File_num)],0,1)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)

    first_guess = [2400, 1e-6, 10, 100e-9]
    
    #Display     
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=40)
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

    #Cleaned - inductive voltage
    ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float), 
            signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*
                            (All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV),
             linewidth = 0.25, alpha = 0.85-(j/50), linestyle = "-", color = "black") #label = "Raw - V(ind)",
   
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", 
            linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 1-(j/50))

    
    #Display crical voltage and current
    #ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
    #        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
    #        linewidth = 3, linestyle = ":", color = "red", alpha = 1, 
    #        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV, Degradaton = " + str(round((100*(1-(popt[0]/Ic_RG))),2)) + "%")
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 25)#, ncols = 4)
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]
    
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + ", IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A]," +  "uV, Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
    
ax.legend(fontsize = 30)#, ncols = 4)

    

#%% Fitting for joint resistance: Ch 1
all_files.sort()
#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 1
Inv_taps = 1
#Starting current
I_start = 300
#End of resisitve Current Value
num_files = len(all_files)

#ch_no = 6
R_ind = 2500
#End of noise Current Value
N_ind = 500

first_guess = [40e-9, 1e-5]
try: 
    ax.cla()
except:
    j = 0
for j in range(0, len(all_files)):
    File_num = F_start + j
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[int(File_num)],0,1)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[int(File_num)],0,1)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]
                             *(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)
 
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))
    fit_y = func_only_R(fit_x, *popt)
    
    fig = plt.figure(1)
    
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])
    fname = fname1[:-4]
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=40)
       
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
    ax.plot(x_data,y_data, linewidth = 0.25, linestyle = "-", color = "black", alpha = 0.95-(j/100))
    ax.plot(fit_x,fit_y,
            label = Tap_name[ch_no] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
            linewidth = 5, linestyle = "--", color = "tab:red")
    
    ax.legend(fontsize = 30)


#%% Graphing Ic against pressure
all_files.sort()
F_start = 0
ch_no = 3
Inv_taps = 1
Mag_f = 1e-7
I_start = 150
num_files = len(all_files)
R_ind = 4499
N_ind = 2500
IvsV = 1
Ic_P = np.zeros([len(all_files),2])

first_guess = [3500, 1e-6, 11, 50e-9]
try: 
    ax.cla()
except: 
    j = 0
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[int(File_num)],0,1)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[int(File_num)],0,1)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    #Display     
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    #fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    fig.suptitle(Tap_name[ch_no] + "- Critical Current vs Pressure (0331-0408)" , fontsize = 40)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_axisbelow(True)
    #ax.set_ylim([3900,4150])
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Pressure [MPa]", fontsize=30)
    ax.set_ylabel("Critical Current [A]", fontsize=30)
    
    #Compression vs IC)
    ax.scatter(All_files[int(File_num)].iloc[0,5],popt[0], s=500, color = "black", marker = "v") #label = fname1

    
    #Save Ic
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]
#ax.scatter(All_files[int(File_num)].iloc[0,5],popt[0], s=500, color = "black", marker = "v", label = Tap_name[ch_no])
Err_std = np.std(Ic_P[:,0], ddof=1) / np.sqrt(np.size(Ic_P[:,0]))

for i in range(0,len(all_files)):
    ax.errorbar(Ic_P[i,1],Ic_P[i,0],Err_std, color = 'black', capsize = 10)

Ic_Pnot0 = np.where(Ic_P[:,1]>10)
Ic_P0 = np.where(Ic_P[:,1]<10)  
Avg_comp = np.mean(Ic_P[Ic_Pnot0[0],0])
Avg_free = np.mean(Ic_P[Ic_P0[0],0])

ax.axhline(y = Avg_free, color = 'tab:red', linestyle = '--', label = "Average Free Ic = " +str(round(Avg_free,2)) + " A -  " + Tap_name[ch_no], lw = 5)
ax.axhline(y = Avg_comp, color = 'tab:blue', linestyle = '--', label = "Average Compressed Ic = " +str(round(Avg_comp,2)) + " A - " + Tap_name[ch_no], lw = 5)
        

#ax.scatter(All_files[int(File_num)].iloc[0,5],popt[0], label = "Ch. no: " + str(ch_no) + " Whole Core", s=500, color = "tab:red")
ax.legend(fontsize = 30)
#%% Graphing Ic against pressure CYCLES

all_files.sort()
F_start = 0
ch_no = 3
Inv_taps = 1
I_start = 150
num_files = len(all_files)
R_ind = 4499
N_ind = 2500
Ic_P = np.zeros([len(all_files),2])

first_guess = [3500, 1e-6, 11, 50e-9]
try: 
    ax.cla()
except: 
    j = 0
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_Noise,All_files[int(File_num)],0,1)
    Avg_at_ResistiveRange_per_tap = Mag_f*Tap_dist[ch_no]*average_in_range(I_indices_R,All_files[int(File_num)],0,1)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*Tap_dist[ch_no]*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    #Display     
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    #fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    fig.suptitle(Tap_name[ch_no] + " - Critical Current vs Total Pressure Cycles (" + fname[0:8] +")" , fontsize = 40)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_axisbelow(True)
    ax.set_ylim([4210,4270])
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Total number of Cycles", fontsize=30)
    ax.set_ylabel("Critical Current [A]", fontsize=30)
    
    #Cycle count vs Ic
    ax.scatter(Cycle_count[j+1],popt[0], s=500, marker = "v", color = "black") #label = fname1
    
    #Save Ic & Pressure
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = j
    
Err_std = np.std(Ic_P[:,0], ddof=1) / np.sqrt(np.size(Ic_P[:,0]))

for i in range(0,len(all_files)):
    ax.errorbar(Cycle_count[i+1],Ic_P[i,0],Err_std, color = 'black', capsize = 10)
    #ax.annotate(Cycle_Pressure[i+1]+"MPa",([Cycle_count[i+1],Ic_P[i,0]]),xytext=(Cycle_count[i+1],Ic_P[i,0]+5),fontsize = 25, 
    #            arrowprops=dict(facecolor='black', arrowstyle="wedge, tail_width = 0.25", alpha=0.1))
    
ax.axhline(y = np.mean(Ic_P[0:4,0]), color = 'tab:red', linestyle = '--', label = "Average Ic before cycling = " +str(round(np.mean(Ic_P[0:4,0]),2)) + " A", lw = 5, alpha = 0.5)
   

#ax.scatter(All_files[int(File_num)].iloc[0,5],popt[0], label = "Ch. no: " + str(ch_no) + " Whole Core", s=500, color = "tab:red")
ax.legend(fontsize = 40)