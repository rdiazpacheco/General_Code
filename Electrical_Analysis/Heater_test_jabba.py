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
import numpy
from scipy.signal import chirp, find_peaks, peak_widths

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

def find_powered_sections(min_height, on_duration,data,columns_per_file):    
    li = []
    for i in range(0,int(len(all_files)/columns_per_file)):
        
        On_indx, pk_hgts = find_peaks(data.iloc[:,columns_per_file*i+1],5)
        li.append(pd.Series(On_indx))
        li[3*i].rename("Off_set " + str(i),inplace = True)
        li.append(pd.Series(On_indx + on_duration))
        li[3*i+1].rename("On_set " + str(i),inplace = True)
        li.append(pd.Series(pk_hgts["peak_heights"]))
        li[3*i+2].rename("Peak_Height " + str(i),inplace = True)
        frame = pd.concat(li, axis=1, ignore_index= False)   
    
    return frame

def find_average_voltages(indices_of_powered_sections, data,columns_per_file):
    li = []
    for j in range(0,int(len(all_files)/columns_per_file)):
        li.append(pd.Series(np.mean(indices_of_powered_sections.iloc[:,3*j+2])))
        li[int(((1/columns_per_file)*len(all_files))*j)].rename("Average_max " + str(j+1),inplace = True)
        lii = np.zeros(len(indices_of_powered_sections.iloc[:,3*j]))
        for i in range(0,len(indices_of_powered_sections.iloc[:,3*j])):
            lii[i] = np.mean(data.iloc[:,2*j+1].iloc[int(indices_of_powered_sections.iloc[:,3*j+1].iloc[i]):int(indices_of_powered_sections.iloc[:,3*j].iloc[i])])
        li.append(pd.Series(lii))
        li[int(((1/columns_per_file)*len(all_files))*j+1)].rename("Average_of_ramp " + str(j+1),inplace = True)

    frame = pd.concat(li, axis=1, ignore_index= False)
    return frame



def Extract_TVI(list_of_files):
    li = []
    for i in range(0,len(list_of_files)):
        fname1 = fname1 = (list_of_files[i].partition('\\')[2])
        fname = fname1[:-6]
        df = pd.read_table(list_of_files[i], delim_whitespace=True)        
        li.append(df.iloc[:,0]) #time stamp
        li.append(df.iloc[:,1]) #voltgae
        li.append(df.iloc[:,7]) # current
        li[3*i].rename("Time: "+ fname, inplace=True)    
        li[3*i+1].rename("Voltage: "+ fname, inplace=True)
        li[3*i+2].rename("Current: "+ fname, inplace=True)
        frame = pd.concat(li, axis=1, ignore_index= False)

    return frame

#%% Extract all data from a day
folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
Vs_1 = Extract_Voltages_only(all_files,1) #extract all data

on_time = 80
OO_indx = find_powered_sections(5,-on_time,Vs_1,2)
V_avgs = find_average_voltages(OO_indx,Vs_1,2)


#%%

test_duration = 3600
for k in range(0,int(0.5*len(all_files))):
    fig, ax = plt.subplots(3,1, figsize = (30,20),constrained_layout = True)
    fname1 = (all_files[2*k+1].partition('\\')[2])
    fname = fname1[:-8]
    Current = fname1[9:-17]
        
    #general set up     
    fig.suptitle(fname1[:-8],fontsize=25)
    for i in range(0,3):
        ax[i].tick_params(axis='x', labelsize=15)
        ax[i].tick_params(axis='y', labelsize=15)
    
    #GRAPH 1
    #Graphing the whole sequence with the ramp points: \
    ax[0].set_title("Whole run", fontsize = 20)
    ax[0].set_ylabel("Voltage [V]", fontsize = 15)
    ax[0].set_xlabel("Time [s]", fontsize = 15)
    ax[0].set_ylim(-1,int(max(Vs_1.iloc[:,2*k+1]))+2)
    
    #Data from repeat 
    time = np.linspace(1,test_duration,len(Vs_1.iloc[:,2*k+1]))
    ax[0].plot(time,Vs_1.iloc[:,2*k+1],linewidth = 1)
    ax[0].plot(OO_indx.iloc[:,3*k]/(len(Vs_1.iloc[:,2*k+1])/test_duration),Vs_1.iloc[:,2*k+1].iloc[OO_indx.iloc[:,3*k]],"x", label = "Peaks of Voltage")
    for i in range(0,len(OO_indx)):
        ax[0].axvline(x = OO_indx.iloc[:,3*k+1].iloc[i]/(len(Vs_1.iloc[:,2*k+1])/test_duration), ymin = 0.05, ymax = 0.95, linestyle = '--', color = 'r',)
        ax[0].axvline(x = OO_indx.iloc[:,3*k].iloc[i]/(len(Vs_1.iloc[:,2*k+1])/test_duration), ymin = 0.05, ymax = 0.95, linestyle = '--', color = 'r',)
    #ax[0].legend(fontsize = 20)
        
    #GRAPH 2
    #Graphing all the traces
    ax[1].set_title("Intra-powered", fontsize = 20)
    ax[1].set_ylabel("Voltage [V]", fontsize = 15)
    ax[1].set_xlabel("Time [s]", fontsize = 15)
    
    for i in range(0,len(OO_indx)):
        ramp_time = np.linspace(0,on_time/10,abs(int(OO_indx.iloc[:,3*k+1].iloc[i]-OO_indx.iloc[:,3*k].iloc[i])))
        ax[1].plot(ramp_time,Vs_1.iloc[:,2*k+1].iloc[int(OO_indx.iloc[:,3*k+1].iloc[i]):int(OO_indx.iloc[:,3*k].iloc[i])])
    
    #GRAPH 3
    #Graph the resistance vs cycles
    ax[2].set_title("Resistance through Cycles", fontsize = 20)
    ax[2].set_ylabel("Resistance [Ohms]", fontsize = 15)
    ax[2].set_xlabel("Cycle Number", fontsize =15)
    ax[2].set_ylim(max(V_avgs.iloc[:,2*k+1]/int(Current))-0.05,max(V_avgs.iloc[:,2*k+1]/int(Current))+0.05)
    
    cycles = np.arange(len(OO_indx.iloc[:,3*k+2]))
    ax[2].scatter(cycles,OO_indx.iloc[:,3*k+2]/int(Current), label = "Max resistance per cycle", s = 60)
    ax[2].axhline(y= np.round(V_avgs.iloc[0,2*k]/int(Current),3), 
                  label = "Avg Highest Resistance: " + str(np.round(V_avgs.iloc[0,2*k]/int(Current),2)) + " Ohms",
                  linestyle  = "--", color = "tab:red", linewidth = 3)
    ax[2].scatter(cycles,V_avgs.iloc[:,2*k+1]/int(Current), label = "Avg resiatance in-cycle", s = 60)
    ax[2].axhline(y= np.round(np.mean(V_avgs.iloc[:,2*k+1]/int(Current)),3), 
                  label = "Avg resistance in-cycle: " + str(np.round(np.mean(V_avgs.iloc[:,2*k+1]/int(Current)),2)) + " Ohms", 
                  linestyle  = "--", color = "tab:orange", linewidth = 3)
    
    ax[2].legend(fontsize = 20)
    
    plt.show()
 #%%
    plt.savefig(fname1[:-8] + ".png", dpi=200)
    
    



    
#%%

folder_path = filedialog.askdirectory()
all_files = glob.glob(os.path.join(folder_path, "*.dat"))

data_strain = Extract_TVI(all_files)
#%%

on_time = 80
On_Off_indexes = find_powered_sections(5,-on_time,data_strain,3)
#V_avgs = find_average_voltages(On_Off_indexes,data_strain,3)


#%% Find peaks
for i in range(0,len(all_files)):
    fname1 = (all_files[i].partition('\\')[2])
    plt.figure(figsize = (15,15))
    peaks, heights = find_peaks(data_strain.iloc[:,3*i+1], height = 1, prominence = 0.5)
    plt.title(fname1[:-8], fontsize = 20)
    
    time = np.arange(len(data_strain.iloc[:,3*i+1]))
    
    plt.plot(time,data_strain.iloc[:,3*i+1])
    plt.scatter(time[peaks],data_strain.iloc[peaks,3*i+1], marker = 'x', color = 'tab:red')
    #plt.scatter(np.arange(len(data_strain.iloc[peaks,1])),data_strain.iloc[peaks,1])
    plt.show()
    
    
#%% comparing Strain vs no strain
li = []
for i in range(0,len(all_files)):
    li.append(all_files[i].partition('\\')[2][:-6])
all_sample_names = pd.DataFrame(li)


#%%




#%%
fig, ax = plt.subplots(3,1, figsize = (20,20))
fname1 = (all_files[0].partition('\\')[2])

fig.suptitle(fname1[:-8],fontsize=25)
ax[0].set_title('Voltage vs Time')
ax[1].set_title('Current vs Time')
ax[2].set_title('Current vs Time')
ax[0].set_ylabel("Volts")
ax[1].set_ylabel("Current")
ax[2].set_ylabel("Volts")
ax[0].set_xlabel("Time")
ax[1].set_xlabel("Time")
ax[2].set_xlabel("Current")
ax[0].plot(data_strain.iloc[:,0],data_strain.iloc[:,1])
ax[1].plot(data_strain.iloc[:,0],data_strain.iloc[:,2])
ax[2].plot(data_strain.iloc[:,2],data_strain.iloc[:,1])
plt.show()

#%%
Vs_1 = data_strain
OO_indx = On_Off_indexes
test_duration = 3600
#for k in range(0,int(0.5*len(all_files))):
for k in range(0,int(2)):
    fig, ax = plt.subplots(3,1, figsize = (30,20),constrained_layout = True)
    fname1 = (all_files[2*k+1].partition('\\')[2])
    fname = fname1[:-8]
    Current = fname1[9:-17]
        
    #general set up     
    fig.suptitle(fname1[:-8],fontsize=25)
    for i in range(0,3):
        ax[i].tick_params(axis='x', labelsize=15)
        ax[i].tick_params(axis='y', labelsize=15)
    
    #GRAPH 1
    #Graphing the whole sequence with the ramp points: \
    ax[0].set_title("Whole run", fontsize = 20)
    ax[0].set_ylabel("Voltage [V]", fontsize = 15)
    ax[0].set_xlabel("Time [s]", fontsize = 15)
    ax[0].set_ylim(-1,int(max(Vs_1.iloc[:,2*k+1]))+2)
    
    #Data from repeat 
    time = np.linspace(1,test_duration,len(Vs_1.iloc[:,2*k+1]))
    ax[0].plot(time,Vs_1.iloc[:,2*k+1],linewidth = 1)
    ax[0].plot(OO_indx.iloc[:,3*k]/(len(Vs_1.iloc[:,2*k+1])/test_duration),Vs_1.iloc[:,2*k+1].iloc[OO_indx.iloc[:,3*k]],"x", label = "Peaks of Voltage")
    for i in range(0,len(OO_indx)):
        ax[0].axvline(x = OO_indx.iloc[:,3*k+1].iloc[i]/(len(Vs_1.iloc[:,2*k+1])/test_duration), ymin = 0.05, ymax = 0.95, linestyle = '--', color = 'r',)
        ax[0].axvline(x = OO_indx.iloc[:,3*k].iloc[i]/(len(Vs_1.iloc[:,2*k+1])/test_duration), ymin = 0.05, ymax = 0.95, linestyle = '--', color = 'r',)
    #ax[0].legend(fontsize = 20)
    