# -*- coding: utf-8 -*-
"""
Jabba Parser Dependencies
Created on Wed Jan 25 14:32:05 2023

@author: rdiazpacheco
"""

# Dependencies

import csv
import time
import json
import copy
import pandas as pd
from tkinter import filedialog
import scipy
#from influxdb import InfluxDBClient
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
from datetime import datetime

#% Functions

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def Extract_voltages_one_file(file_index,j,filename,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue,load_cell):   
    li = []
    df = pd.read_csv(filename, header=header_no, skiprows = range(7,24))
    fname1 = (file_index[j].partition('\\')[2])
    fname = fname1[:-4]
    try:
        I_val = df.loc[:,I_col_name];
        
    except: 
        header_no = 4
        df = pd.read_csv(filename, header=header_no, skiprows = range(7,24))
        fname1 = (file_index[j].partition('\\')[2])
        fname = fname1[:-4]
        I_val = df.loc[:,I_col_name];
        
    li.append(I_val)
    li[0].rename("I_"+fname,inplace=True)

    for i in range(1,(1+no_ch)):
        try: 
            VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
            li.append(VTap)
            li[i].rename("V"+str(i),inplace=True)
        except: 
            li.append(pd.Series(0))
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

def offset_voltage_perCh(data,IvsV,total_number_ch):
    data = data.replace(r'^\s*$', np.nan, regex=True)
    I_start = np.where(data.iloc[:,0].astype(float)>0)
    I_start_ind = min(I_start[0])
    Offset_values = []
    for i in range(0,total_number_ch):
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

def average_in_range(indices,data,total_number_ch,filt_before,IvsV):  
    Avg_noise1 = []
    if filt_before == 1:    
        for i in range(0,total_number_ch):
            Avg_noise1.append(np.mean(
                butter_lowpass_filter(cutoff, fs, order,data.iloc[int(IvsV*indices[0]):int(IvsV*indices[1]),i+1].astype(float))))           
    else:
        for i in range(0,total_number_ch):
            Avg_noise1.append(np.mean(data.iloc[int(IvsV*indices[0]):int(IvsV*indices[1]),i+1].astype(float)))
    Avg_noise = pd.Series(Avg_noise1)                  
    return Avg_noise

def average_value_at_step(data,step_currents,decay_time,IvsV,tot_ch,Imax):
    #step_currents is an array with the step values
    data = data.replace(r'^\s*$', np.nan, regex=True)
    average_Vstep_per_tap = []
    aaa= []
    for i in range(0,len(step_currents)):
        I_nA = []
        I_nB = []   
        I_ABs = np.zeros(2)
        I_nAp = np.where(data.iloc[:,0].astype(float) > step_currents[i]-1)
        I_nA.append(I_nAp[0])   
        I_ABs[0] = int(min(I_nA[0]))
        if step_currents[i] == Imax:
            I_nBp = np.where(data.iloc[:,0].astype(float) == step_currents[i])
            I_nB.append(I_nBp[0])
            I_ABs[1] = int(max(I_nB[0]))
        else: 
            I_nBp = np.where(data.iloc[:,0].astype(float) > step_currents[i]+1)
            I_nB.append(I_nBp[0])
            I_ABs[1] = int(min(I_nB[0]))
        range_indices = I_ABs
        average_in_step = []
        
        for j in range(1,tot_ch+1):
            average_in_step.append(np.mean(data.iloc[int((IvsV*range_indices[0])+decay_time):int(IvsV*range_indices[1]),j]))
        aaa.append(pd.Series(average_in_step))
        average_Vstep_per_tap = pd.concat(aaa, axis=1, ignore_index= False) 
        
    return average_Vstep_per_tap


def func(x, Vc, Ic, V_floor, n):
    return Vc*(x/Ic)**n + V_floor

def func_w_R(x, Vc, Ic, V_floor, n,R):
    return Vc*(x/Ic)**n + V_floor + x*R

def func_only_R(x,R,OffSet):
    return x*R+OffSet

def assign_names(n1, n2, n3, n4):
    vtap_names = []
    vtap_names.append(str(n1))
    vtap_names.append(str(n2))
    vtap_names.append(str(n3))
    vtap_names.append(str(n4))
    return vtap_names
