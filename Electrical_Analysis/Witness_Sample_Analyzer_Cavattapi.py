# -*- coding: utf-8 -*-
"""
Witness Samples
Created on Wed Dec 21 12:39:30 2022

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

def Extract_voltages_one_file2(filename,header_no,I_col_name,V_cols_name):
    li = []
    df = pd.read_csv(filename, header=header_no, skiprows = range(7,24))
    fname1 = (all_files[j].partition('\\')[2])
    fname = fname1[-7:-4]
    I_val = df.loc[:,I_col_name];
    li.append(I_val)
    li[0].rename("I_"+fname,inplace=True)
    for i in range(1,(1+len(V_cols_name))):
        VTap = df.loc[:,str(V_cols_name[i]) + " [uV]"]
        li.append(VTap)        
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
#%% Data import - pre
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
all_files.sort()
num_files = len(all_files)
indx1 = []
for j in range(0,num_files):
    indx1.append(int(all_files[j].partition('\\')[2][-7:-4]))    
all_files1 = pd.DataFrame(all_files,index = indx1)

# Data import - post
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
all_files.sort()
num_files = len(all_files)
indx1 = []
for j in range(0,num_files):
    indx1.append(int(all_files[j].partition('\\')[2][-7:-4]))    
all_files2 = pd.DataFrame(all_files,index = indx1)

#%% Indicate files of interest 

Cavatappi_Pre0 = {
    "Pb1":[41,42,43],
    "Pb2":[57,58,59],
    "Pb3":[64,65],
    "Pb4":[71,72,73],
    "Pb5":[76,77,78],
    "Ni1":[81,82,83],
    "Ni2":[88,89,90],
    "Ni3":[92,93,94],
    "Ni4":[96,97,98],
    "Ni5":[100,101,102]   
    }
Cavatappi_Pre = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in Cavatappi_Pre0.items() ]))


Cavatappi_Post0 = {
    "Pb1":[15,16],
    "Pb2":[28,29,30],
    "Pb3":[34],
    "Pb4":[38,39,40],
    "Pb5":[48,49,50],
    "Ni1":[54,55,56],
    "Ni2":[60,61,62],
    "Ni3":[65,66,67],
    "Ni4":[67,70,71],
    "Ni5":[73,74,75]    
    }
Cavatappi_Post = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in Cavatappi_Post0.items() ]))


tap_names = {
    1:"5cm",
    2:"10cm"
    }
tap_names2 = {
    1:"5_cm",
    2:"10_cm"
    }
tot_ch = 2
  
#% Extract voltages: One dictionary by sample
All_files1 = {}
All_files2 = {}
Sample_name = "Pb1"
#For Pre
for k in range(0,np.count_nonzero(~np.isnan(Cavatappi_Pre.loc[:,Sample_name]))):  
    one_filename = all_files1.loc[int(Cavatappi_Pre[Sample_name][k]),0]
    All_files1[int(one_filename[-7:-4])] = Extract_voltages_one_file2(one_filename,5,"Current [A]",tap_names)
    
#For Post
for k in range(0,np.count_nonzero(~np.isnan(Cavatappi_Post.loc[:,Sample_name]))):  
    one_filename = all_files2.loc[int(Cavatappi_Post[Sample_name][k]),0]
    All_files2[int(one_filename[-7:-4])] = Extract_voltages_one_file2(one_filename,5,"Current [A]",tap_names2)


#% Fitting without steps all in the file
#""" 20221116
Tap_dist = {
    1:5,
    2:10
    }
#"""

#""" 20221116
Tap_name = {
    1:"5 cm",
    2:"10 cm"
    }
#"""'

I_start = 5
#End of resisitve Current Value
num_files = np.count_nonzero(~np.isnan(Cavatappi_Pre.loc[:,Sample_name]))
ch_no = 1
R_ind = 100
#End of noise Current Value
N_ind = 159

#"""
for k in range(0,np.count_nonzero(~np.isnan(Cavatappi_Pre.loc[:,Sample_name]))): 
    Current_indices, Imax = find_start_end_ramp_onefile(All_files1[int(Cavatappi_Pre[Sample_name][k])],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files1[int(Cavatappi_Pre[Sample_name][k])],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files1[int(Cavatappi_Pre[Sample_name][k])],I_start, Imax-1)
    Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files1[int(Cavatappi_Pre[Sample_name][0])],0)
    
    x_data = All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e0))):
        return Vc*(x/Ic)**n + V_floor + x*R
    first_guess = [120, 1, 10, 1]
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,Imax,len(All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)
    
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    #fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    #fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Before vs After Solder, Witness Sample: " + Sample_name,fontsize=25)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    """
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.000025))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)
    

    #Raw - inductive voltage
    ax.plot(All_files1[int(Cavatappi_Pre[Sample_name][0])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            All_files1[int(Cavatappi_Pre[Sample_name][0])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1],
            label = "Raw - V(ind)", linewidth = 1, linestyle = "-", color = "tab:blue", alpha = 0.9-(1/10))
    
    
    #Fitted Function
    ax.plot(fit_x, fit_y,label = "BEFORE SOLDER Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1))
            + " [A], R=" + str(round(popt[3]*(1e2),1)) + " n\u03A9 ", linewidth = 7, linestyle = "-.", color ="tab:green", alpha = 0.9-(1/10))
    
    
    #Display crical voltage and current
    ax.plot(All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            np.full_like(All_files1[int(Cavatappi_Pre[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e0),
            linewidth = 3, linestyle = ":", color = "blue", alpha = 1, label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
    ax.axvline(x = round(popt[0],1), color = 'green', linewidth = 2,linestyle='-.')
    
#"""

#POST SOLDER
ch_no = 1
N_ind = Imax-1

for k in range(0,np.count_nonzero(~np.isnan(Cavatappi_Post.loc[:,Sample_name]))): 
    Current_indices, Imax = find_start_end_ramp_onefile(All_files2[int(Cavatappi_Post[Sample_name][k])],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files2[int(Cavatappi_Post[Sample_name][k])],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files2[int(Cavatappi_Post[Sample_name][k])],I_start, Imax-1)
    Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files2[int(Cavatappi_Post[Sample_name][0])],0)
    
    x_data = All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0]
    y_data = All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e0))):
        return Vc*(x/Ic)**n + V_floor + x*R
    first_guess = [100, 0, 21, 1]
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,Imax,len(All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no]))
    fit_x2 = np.linspace(0,Imax,len(All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)
    
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    #fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    #fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    #fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    """
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.000025))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)
    
    #Raw - inductive voltage
    ax.plot(All_files2[int(Cavatappi_Post[Sample_name][0])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            All_files2[int(Cavatappi_Post[Sample_name][0])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1],
            label = "Raw - V(ind)", linewidth = 1, linestyle = "-", color = "tab:orange",alpha = 0.9-(1/10))
    
    
    #Fitted Function
    ax.plot(fit_x, fit_y,label = "AFTER SOLDER Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e2),1)) + " n\u03A9 ", linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 0.9-(1/10))
    
    
    #Display crical voltage and current
    ax.plot(All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],
            np.full_like(All_files2[int(Cavatappi_Post[Sample_name][k])].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e0),
            linewidth = 3, linestyle = ":", color = "red", alpha = 1, label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
ax.legend(fontsize = 15) 