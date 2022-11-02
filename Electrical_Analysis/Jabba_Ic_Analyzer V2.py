# -*- coding: utf-8 -*-
"""
Jabba Ic Analyzer - more robust parser
Created on Wed Oct  5 16:45:31 2022

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

#% Functions

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def Extract_Voltages(all_files,header_no,I_col_name,V_cols_name,rvalORvalue):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):        
        df = pd.read_csv(all_files[j], header=header_no, skiprows = range(7,24))
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        I_val = df.loc[:,I_col_name];
        li.append(I_val)
        li[(5*j)].rename("I_"+fname,inplace=True)
        for i in range(1,5):
            VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
            li.append(VTap)
            li[(i)+(5*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)            
        frame = pd.concat(li, axis=1, ignore_index= False)        
    return frame

def Extract_Voltages(all_files,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):        
        df = pd.read_csv(all_files[j], header=header_no, skiprows = range(7,24))
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        I_val = df.loc[:,I_col_name];
        li.append(I_val)
        li[((1+no_ch)*j)].rename("I_"+fname,inplace=True)
        for i in range(1,(1+no_ch)):
            VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
            li.append(VTap)
            li[(i)+((1+no_ch)*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)            
        frame = pd.concat(li, axis=1, ignore_index= False)        
    return frame

def get_average_of_tap(files_name,data,Current_indices,tap_of_interest):
    li = [] 
    li.append(pd.Series(data.iloc[:,0].iloc[int(Current_indices[0,0]):int(Current_indices[0,1])].values))

    li1 = []
    for j in range(0,len(files_name)):
        li1.append(data.iloc[:,(1+no_ch)*j+(tap_of_interest)].iloc[int(Current_indices[j,0]):int(Current_indices[j,1])])
        li2 = np.mean(li1, axis = 0)
    li.append(pd.Series(li2))
       
    avg_V_tap = pd.concat(li, axis=1, ignore_index= False)    
    
    return avg_V_tap

def resistive_range (files_name,data,Istart,Iend):
    I_nA = []
    I_nB = []   

    I_ABs = np.zeros([len(files_name),2])
    for i in range(0,len(files_name)):
        I_nAp = np.where(data.iloc[:,(1+no_ch)*i] > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[i,0] = min(I_nA[i])
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(InVs_1.iloc[:,(1+no_ch)*i]):  
            I_nB.append(min(I_nA[i]) + 500)
            I_ABs[i,1] = min(I_nA[i]) + 500
        else:
            I_nBp = np.where(data.iloc[:,(1+no_ch)*i] > Iend)
            I_nB.append(I_nBp[0])    
            I_ABs[i,1] = min(I_nB[i])          
    return I_ABs

def noise_range(files_name,data,Istart,Iend):
    I_nA = []
    I_nB = []   

    I_ABs = np.zeros([len(files_name),2])
    for i in range(0,len(files_name)):
        I_nAp = np.where(data.iloc[:,(1+no_ch)*i] > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[i,0] = min(I_nA[i])
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(InVs_1.iloc[:,(1+no_ch)*i]):  
            I_nB.append(min(I_nA[i]) + 500)
            I_ABs[i,1] = min(I_nA[i]) + 500
        else:
            I_nBp = np.where(data.iloc[:,(1+no_ch)*i] > Iend)
            I_nB.append(I_nBp[0])    
            I_ABs[i,1] = min(I_nB[i])          
    return I_ABs

def find_noise_per_tap(Noice_indices,data,filt_before):
  
    Avg_noise1 = []
    if filt_before == 1:    
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise1.append(np.mean(
                    butter_lowpass_filter(cutoff, fs, order,data.iloc[int(Noice_indices[j,0]):int(Noice_indices[j,1]),j*(1+no_ch)+i+1])))   
    else:
        for j in range(0,len(all_files)):
            for i in range(0,4):
                Avg_noise1.append(np.mean(data.iloc[int(Noice_indices[j,0]):int(Noice_indices[j,1]),j*(1+no_ch)+i+1]))
    Avg_noise = pd.Series(Avg_noise1)         
    
    avg_noise_p_tap = []           
    for i in range(0,len(all_files)):
        avg_noise_p_tap.append(np.mean(Avg_noise[[i,i+4,i+8,i+12]]))
      
    return Avg_noise, avg_noise_p_tap

def find_start_end_ramp(data,Istart):
    I_ind = []
    max_Iall = []
    I_str_stp = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nparray = np.where(data.iloc[:,5*i]>Istart)
        I_ind.append(I_nparray[0])
        I_max = max(data.iloc[:,5*i])
        max_Ip = np.where(data.iloc[:,5*i] == max(data.iloc[:,5*i]))
        max_Iall.append(max_Ip[0])
        I_str_stp[i,0] = min(I_ind[i])
        I_str_stp[i,1] = min(max_Iall[i])        
    return I_str_stp, I_max

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

def func_w_R(x, Vc, Ic, V_floor, n,R):
    return Vc*(x/Ic)**n + V_floor + Vc*x*R

def func_only_R(x, Vc,R):
    return Vc*x*R

def assign_names(n1, n2, n3, n4):
    vtap_names = []
    vtap_names.append(str(n1))
    vtap_names.append(str(n2))
    vtap_names.append(str(n3))
    vtap_names.append(str(n4))
    return vtap_names


#%%
# Select a folder and then extract the data into one Data Frame
# I_col_name refers to the name of the Current Column in the CSV file
#V_cols_name refers to the name of the Voltage taps. Could be ChStat[1].Value or VTapFilter[1].rVal - we should change this on Jabba
no_ch = 8
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
InVs_1 = Extract_Voltages(all_files,8,4,"CmdAmps.Value","ChStat","Value")

#%%
#butter filter specs
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 2      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

#%%
R_ind = 1700

#Select noice range, give it a intial and final current where the noise should be
Current_indices, Imax = find_start_end_ramp(InVs_1,200)
Noice_ind = noise_range(all_files,InVs_1,200, 400)
Resitive_ind = noise_range(all_files,InVs_1,200, R_ind)

#%% Select the values where the noise should be. This may or may not be used later
Avg_noises, n_p_t_avg = find_noise_per_tap(Noice_ind,InVs_1,0)

#%%
Ic_vals, Ic_avg = find_Ic_values(1,1,1,1,Avg_noises)





#%% Display one tap a time with average values and fittings
sample_name = "Queen Snake, Sample 2: QS_A1-QS_A2: "
tap_of_interest = 1
n_guess = 5

vtns = assign_names("Lead QS_A1", "Cable: QS_A2", "PIT-V Joint", "Lead QS_A2") #for 20220819
 

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
                       InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                       label="Voltage Tap", linewidth = 3,c = c)
            #InVs_1["VTap"+str(i+1+2*j)+"_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*k)+(2*j)+i]
            c = next(color) 
            c = next(color) 
            IV_1, = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                        butter_lowpass_filter(cutoff, fs, order,InVs_1["VTap"+str(i+1+2*j)+"_"+fname]
                                              .iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]),
                        label="2Hz LPF", linewidth = 3,c = c)
            #.iloc[int(Current_indices[k,0]):int(Current_indices[k,1])]-Avg_noises[(4*k)+(2*j)+i]),
            
            # 1 uV/cm line
            OneuV = ax[i,j].plot(InVs_1["I_"+fname].iloc[int(Current_indices[k,0]):int(Current_indices[k,1])],
                                 SC_parameter, label="1uV/cm", color='red', linewidth=4, linestyle=':')
            #fit functions 
            first_guessRO = [1, 0.001] # Linear fit Only
            first_guessFULL = [1, 2000, 0.1, 10, 0.01 ] # full fit
            popt, pcov = curve_fit(func_only_R,InVs_1.iloc[int(Resitive_ind[0,0]):int(Resitive_ind[0,1]),0], 
                                   InVs_1.iloc[int(Resitive_ind[0,0]):int(Resitive_ind[0,1]),i+j+1], 
                                   p0=first_guessRO)
            popt2, pcov2 = curve_fit(func_w_R,InVs_1.iloc[int(Current_indices[0,0]):int(Current_indices[0,1]),0], 
                                   InVs_1.iloc[int(Current_indices[0,0]):int(Current_indices[0,1]),i+j+1], 
                                   p0=first_guessFULL)
            fit_x = np.linspace(200,R_ind,3000)
            fit_x2 = np.linspace(200,Imax,3000)
            fit_y = func_only_R(fit_x, *popt)
            fit_y2 = func_w_R(fit_x2, *popt2)
            
            
            Fit_lineRO = ax[i,j].plot(fit_x,fit_y)
            Fit_lineR = ax[i,j].plot(fit_x2,fit_y2)

            
    
    for j in range(0,2):
        for i in range(0,2):           
            leg = ax[i,j].legend(fontsize=15)




#%%
# This will focus on one tap and do all the fittings and find Ic and N values (the fitting)
# 1. Select tap number (1,2,3,4)
# Select a guess for the n value, usually something between 5-10. You can reiterate to get a better fit. If it's too far (over or under, the fit will be obviosly bad)
# give the sample a name "queen snake"
# and that's it. It will average all the runs of the same kind of the same filter and do all the analysis.
vtns = assign_names("Lead QS_A1", "Cable: QS_A2", "PITV Joint", "Lead QS_A2") #for 20220819
#Select noice range
#sample_name = "QS 2, Ch 4 "
sample_name = "PIT VIPER Joint QS_A1-QS_A2 "
#sample_name = "Overall QS_A1-QS_A3 "
#sample_name = "QS_A2 Cable "
#sample_name = "QS2_A2 Lead "
tap_of_interest = 6
tap_dist = 25
n_guess = 8

#butter low pass filter is used here
T = 5.0         # Sample Period
fs = 500.0       # sample rate, Hz
cutoff = 1     # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples

#Ic_vals, Ic_avg = find_Ic_values(1,1,1,1,Avg_noises)

#get average run
Avg_InV_tap = get_average_of_tap(all_files, InVs_1, Current_indices, tap_of_interest)

def func_only_R(x,R):
    return x*R
#fit functions 
first_guessRO = [0.0001] # Linear fit Only

#popt, pcov = curve_fit(func_only_R,Avg_InV_tap.iloc[int(Resitive_ind[int(tap_of_interest-1),0]):int(Resitive_ind[int(tap_of_interest-1),1]),0], 
#                       tap_dist*Avg_InV_tap.iloc[int(Resitive_ind[int(tap_of_interest-1),0]):int(Resitive_ind[int(tap_of_interest-1),1]),1], 
#                       p0=first_guessRO)
#fit_x = np.linspace(200,R_ind,3000)
#fit_y = (1/tap_dist)*func_only_R(fit_x, *popt)

first_guessFULL = [1, 1900, 0.0001, 17] # full fit
#first_guessFULL = [1, 50, 1, 20, 0.001] # full fi


popt2, pcov2 = curve_fit(func,Avg_InV_tap.iloc[int(Current_indices[int(tap_of_interest-1),0]):int(Current_indices[int(tap_of_interest-1),1]),0]-n_p_t_avg[int(tap_of_interest)-1], 
                       Avg_InV_tap.iloc[int(Current_indices[int(tap_of_interest-1),0]):int(Current_indices[int(tap_of_interest-1),1]),1]+n_p_t_avg[int(tap_of_interest)-1], 
                       p0=first_guessFULL)
#popt2, pcov2 = curve_fit(func_w_R,Avg_InV_tap.iloc[int(Current_indices[int(tap_of_interest-1),0]):int(Current_indices[int(tap_of_interest-1),1]),0], 
#                      -1*Avg_InV_tap.iloc[int(Current_indices[int(tap_of_interest-1),0]):int(Current_indices[int(tap_of_interest-1),1]),1]+n_p_t_avg[int(tap_of_interest)-1], 
#                       p0=first_guessFULL)

fit_x2 = np.linspace(200,Imax,3000)
#fit_y2 = func_w_R(fit_x2, *popt2)
fit_y2 = func(fit_x2, *popt2)
SC_parameter = np.zeros(len(Avg_InV_tap))
SC_parameter[:] = 1                     

plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()

i=0
fname1 = (all_files[i].partition('\\')[2])
fname = fname1[:-4]
fig.suptitle("Ch. " + str(tap_of_interest) + " - " + sample_name + fname[:-4],fontsize=25)

ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
#plt.figure(figsize = (30,18))
ax.set_ylim((-5,50))
ax.set_xlim((50,int(Imax)))
ax.set_xlabel("Current [A]", fontsize=15)
ax.set_ylabel("Electric Field [uV/cm]", fontsize=15)

ax.plot(Avg_InV_tap.iloc[:,0],-0.1*butter_lowpass_filter(cutoff, fs, order,Avg_InV_tap.iloc[:,1]), 
         label = "Average Raw", linewidth = 2, linestyle = "-", color = "tab:green")
ax.plot(fit_x, fit_y,label = "Lin fit, R = " + str(round(1000*popt[0],4)) + " nOhm", linewidth = 5, linestyle = "-.", color = "red", alpha = 0.7)
#ax.plot(fit_x2, fit_y2,label = "Exp fit, Ic = " + str(round(popt2[1],2)) + " A," + " n = " +str(round(popt2[3],2)), linewidth = 4, linestyle = "-", color = "tab:blue", alpha = 0.7)

#-n_p_t_avg[tap_of_interest-1]
for i in range(0,len(all_files)):
    fname1 = (all_files[i].partition('\\')[2])
    fname = fname1[:-4]
    ax.scatter(InVs_1.iloc[int(Current_indices[i,0]):int(Current_indices[i,1]),5*i],
                           -0.1*InVs_1.iloc[int(Current_indices[i,0]):int(Current_indices[i,1]),5*i+tap_of_interest],
                           label = "Raw " + fname,marker= 10, s= 6, alpha=0.2)
ax.plot(Avg_InV_tap.iloc[:,0],SC_parameter, label="1uV/cm", color='red', linewidth=2, linestyle=':')
ax.axvline(x = round(popt2[1],1), color = 'red', linewidth = 2,linestyle=':')
ax.legend(fontsize = 20)
fig.show()
filename4 = sample_name + fname
#plt.savefig("Lead QS1_A3 Lead "+ fname,dpi = 200)
#%%

plt.savefig(filename4,dpi = 200)



#%%
#def Extract_Voltages2(all_files,header_no,I_col_name,V_cols_name,rvalORvalue):
header_no = 6
I_col_name = "CmdAmps.Value"
V_cols_name = "VTapFilter"
rvalORvalue = "rVal"

li = []
num_files = len(all_files)
for j in range(0,num_files):        
    df = pd.read_csv(all_files[j], header=header_no, skiprows = range(7,24))
    fname1 = (all_files[j].partition('\\')[2])
    fname = fname1[:-4]
    I_val = df.loc[:,I_col_name];
    li.append(I_val)
    li[(5*j)].rename("I_"+fname,inplace=True)
    for i in range(1,5):
        VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
        li.append(VTap)
        li[(i)+(5*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)            
    frame = pd.concat(li, axis=1, ignore_index= False)    
        
#    return frame

