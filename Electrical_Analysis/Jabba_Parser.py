# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 08:32:25 2022
Jabba CSV Output Parser

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

def find_noise_per_tap(I_ABs,filt_before):
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
    
    Avg_noise = []
    if filt_before == 1:    
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise.append(np.mean(butter_lowpass_filter(cutoff, fs, order,InVs_1.iloc[int(I_ABs[j,0]):int(I_ABs[j,1]),j*5+i+1])))   
    else:
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise.append(np.mean(InVs_1.iloc[int(I_ABs[j,0]):int(I_ABs[j,1]),j*5+i+1]))
    return Avg_noise

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def find_Ic_values(a,b,c,d,noise_sub,avg_npt):
    #butter low pass filter is used here
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
    
    if noise_sub == 1:
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
                    
    else:    
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

def find_n_val(a,b,c,d,avg_npt):
    T = 5.0         # Sample Period
    fs = 500.0       # sample rate, Hz
    cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples
    log_InV = []
    for j in range(0,len(all_files)):
        log_InV1 = pd.DataFrame(np.log(InVs_1.iloc[:,5*j]))
        log_InV.append(log_InV1)
        #log_InV[(5*j)].rename("logI_"+fname,inplace=True)
        for i in range(0,4):
            log_InV2 = np.log(abs(butter_lowpass_filter(cutoff, fs, order, InVs_1.iloc[:,(5*j)+i+1].iloc[int(Current_indices[j,0]):int(Current_indices[j,1])]-avg_npt[(4*j)+i])))
            log_InV3 = pd.DataFrame(log_InV2)
            log_InV.append(log_InV3)
            #the vtaps after the filter cannot conserve their native index for now    
        log_InVf = pd.concat(log_InV, axis=1, ignore_index= False) 

    selected_graphs = (a,b,c,d)
    slopes =  np.zeros([4*len(all_files),1])
    for j in range(0,len(all_files)):        
        for i in range(0,4):
            if selected_graphs[i] == 1:
                slope,intercept,rval,pval,stderrrenamestderr = linregress(log_InVf.iloc[:,5*j].iloc[int(Current_indices[j,0]):int(Current_indices[j,1])],
                                       log_InVf.iloc[0:int(Current_indices[j,1])-int(Current_indices[j,0]),((5*j)+i+1)])
                slopes[4*j+i] = slope
            else:
                slopes[4*j+i] = 0
        
    return slopes


def assign_names(n1, n2, n3, n4):
    vtap_names = []
    vtap_names.append(str(n1))
    vtap_names.append(str(n2))
    vtap_names.append(str(n3))
    vtap_names.append(str(n4))
    return vtap_names


def find_fit_func(xdata,ydata,I_exp_start,a,b,c,d):
    selected_graphs = (a,b,c,d)
    exp_coeff = np.zeros((len(selected_graphs),2))
    I_exp_ind = np.where(xdata == I_exp_start)
    
    if I_exp_ind[0][0] > 0:
        for t in range(0,len(selected_graphs)):
            if selected_graphs[t] == 1:
                x = InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]
                y = pd.Series(abs(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])])-Avg_noises[(4*j)+i]))
                a_ex, b_ex = np.polyfit(x[I_exp_start:len(x)], np.log(y[I_exp_start:len(x)]), 1, w=np.sqrt(y[I_exp_start:len(x)]))
                exp_coeff[t,0] = a_ex
                exp_coeff[t,1] = b_ex                
            else:
                exp_coeff[t,0] = 0
                exp_coeff[t,1] = 0                
    else: 
        exp_coeff[:,0] = 0
        exp_coeff[:,1] = 0
    return exp_coeff

def func(x, a, b, c):
    return np.exp(b) * np.exp(a * x) + c
#%% Extract all data from a day
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
InVs_1 = Extract_Voltages(all_files)

# Cleaning Data: subtracting linear noise
Noice_ind = noise_range(200, 800)
Avg_noises = find_noise_per_tap(Noice_ind,1)
#Avg_noises_f = find_noise_per_tap_f(Noice_ind)

# Find Ic and n value
Current_indices = find_start_end_ramp(200)
Ic_values = find_Ic_values(1, 1, 0, 1, 1,Avg_noises)

slopes = find_n_val(1,1,0,1,Avg_noises)




#%%
#20220816_003 abd 019 are not working. Upon opening and saving the csv on excel. Check on text editor next


#%%

#vtns = assign_names("Lead 1", "Lead 2", "PIT-V Joint", "Overall") #for 20220816
#vtns = assign_names("Lead 1", "Lead 2", "PIT-V Joint", "Cable: A1") #for 20220817
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
                       InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*j)+i],
                       label="Voltage Tap", linewidth = 3,c = c)
            c = next(color) 
            c = next(color) 
            IV_1, = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                        butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                                              .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*j)+i]),
                        label="2Hz LPF, IC = "+ str(Ic_values[4*k+2*j+i,1]) + ", n = " + str(slopes[4*k+2*j+i,0])[0:3], linewidth = 3,c = c)
            
            OneuV = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                                 SC_parameter, label="1uV/cm", color='red', linewidth=4, linestyle=':')
    
    for j in range(0,2):
        for i in range(0,2):           
            leg = ax[i,j].legend(fontsize=15)
  


#%% checking different fits
Noice_ind = noise_range(200, 900)
Avg_noises1 = find_noise_per_tap(Noice_ind,0)
Avg_noises2 = find_noise_per_tap(Noice_ind,1)

Current_indices = find_start_end_ramp(200)
Ic_values = find_Ic_values(1, 1, 0, 1, 1,Avg_noises)


j= 0
i=0
k= 5 
fname1 = (all_files[k].partition('\\')[2])
fname = fname1[:-4]

   #butter low pass filter
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

    
test_I_val = InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]


test_V_val = InVs_1["VTap"+str(i+1)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*j)+i]

test_V_val_clean = pd.Series(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]))
test_V_val_clean_minus_noice = pd.Series(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])])-Avg_noises1[(4*j)+i])
test_V_val_clean_minus_noice_f = pd.Series(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises2[(4*j)+i]))
test_0_val = np.zeros((len(test_I_val),1))
test_1_val = np.zeros((len(test_I_val),1))
test_1_val[:] = 1
#%
plt.figure(figsize=(30, 20), dpi=100)
plt.plot(test_I_val,test_V_val,label = "raw")
plt.plot(test_I_val,test_V_val_clean,label = "LPF")
plt.plot(test_I_val,test_V_val_clean_minus_noice,label = "LPF minus noise")
plt.plot(test_I_val,test_V_val_clean_minus_noice_f,label = "LPF minus noise_f")
plt.plot(test_I_val,test_0_val,label = "0 V",color='red', linewidth=4, linestyle=':')
plt.plot(test_I_val,test_1_val,label = "1 V",color='red', linewidth=4, linestyle=':')
plt.title(fname,fontsize=25)


#test_V_val_clean_minus_noice_f
plt.legend(fontsize = 20)
#%% Find exponetial fit
from scipy.optimize import curve_fit
j= 0
i=0
k= 5 
fname1 = (all_files[k].partition('\\')[2])
fname = fname1[:-4]

   #butter low pass filter
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

Noice_ind = noise_range(200, 900)
Avg_noises = find_noise_per_tap(Noice_ind,1)

x = InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]
y = pd.Series(22*abs(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])])-Avg_noises[(4*j)+i]))

#aex2, bex2 = np.polyfit(x, np.log(y), 1, w=np.sqrt(y))
aex2, bex2 = np.polyfit(x[3000:len(x)], np.log(y[3000:len(x)]), 1, w=np.sqrt(y[3000:len(x)]))

def find_fit_func(xdata,ydata,I_exp_start):
    x = InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]
    y = pd.Series(abs(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                          .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])])-Avg_noises[(4*j)+i]))
    a_ex, b_ex = np.polyfit(x[I_exp_start:len(x)], np.log(y[I_exp_start:len(x)]), 1, w=np.sqrt(y[I_exp_start:len(x)]))

    return a_ex, b_ex

def func(x, a, b, c):
    return np.exp(b) * np.exp(a * x) + c


test_I_val = InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]
test_V_val = 22*InVs_1["VTap"+str(i+1)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*j)+i]

test_V_val_clean = 22*pd.Series(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]))
test_V_val_clean_minus_noice = 22*pd.Series(butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                      .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])])-Avg_noises[(4*j)+i])
test_0_val = np.zeros((len(test_I_val),1))
test_1_val = np.zeros((len(test_I_val),1))
test_1_val[:] = 1
#%
plt.figure(figsize=(30, 20), dpi=100)
plt.plot(test_I_val,test_V_val,label = "raw")
plt.plot(test_I_val,test_V_val_clean,label = "LPF", linewidth = 4)
plt.plot(test_I_val,test_V_val_clean_minus_noice,label = "LPF minus noise", linewidth = 4)

#plt.plot(x,func(x,aex1,bex1,0),label="fitted w=0", color = "green", linewidth = 4)
plt.plot(x,func(x,aex2,bex2,0),label="fitted, n = " + "{:.3f}".format(aex2), color = "red", linewidth = 4)
        
plt.plot(test_I_val,test_0_val,label = "0 V",color='red', linewidth=4, linestyle=':')
plt.plot(test_I_val,test_1_val,label = "1 V",color='red', linewidth=4, linestyle=':')
plt.title(fname,fontsize=25)


#test_V_val_clean_minus_noice_f
plt.legend(fontsize = 20)




