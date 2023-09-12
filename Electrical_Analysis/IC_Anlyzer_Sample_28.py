# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:09:46 2023
Analyser for IV data from Graphana

@author: rdiazpacheco
"""

#%%
#% Data import
import os
import csv
import matplotlib.pyplot as plt

os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
folder_path = filedialog.askdirectory()
#folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))

#%%
filename = all_files[0]
all_data1 = pd.read_csv(filename, header=0)#, skiprows = range(7,24))

#Name dictionary 
channels = {
    0:"time",
    1:"V1_V7",
    2:"V2_V4",
    3:"V2_V5",
    4:"V2_V6",
    5:"V3_V5",
    6:"V3_V6",
    7:"V4_V6",
    8:"isrc",
    9:"LN2 Bath Bulk Temp 1",
    10:"LN2 Bath Bulk Temp 2",
    11:"ref_timestamp"
    }

#Tap distance in cm
tap_dist = {
    1:60,
    2:24.5,
    3:37.5,
    4:50.5,
    5:15.0,
    6:38.0,
    7:26.0,
    }
#%% Find where the current is greater than 0

""" 
This should apply to one or many runs in one fil
Make sure to make "shot_section_start" = 0 !!!
"""
shot_section_start = 8000 #select where to start looking at data. For single shots, this would be 0
all_data = all_data1.iloc[shot_section_start:-1,:]

#current_on_indexes = find(all_data.loc[:,"isrc"])
v_tap = 6
current_on_matrix_1 = all_data[all_data.loc[:,"isrc"]*1000 > 100] #convert kA to A and cut anything beloow 500
current_on_matrix = current_on_matrix_1[current_on_matrix_1.loc[:,channels[v_tap]]>0] # we get only the data [voltage] - this way we eliminate some of the ramp downs. 

#%% if running only on
all_runs = {}
all_runs[0] = current_on_matrix

#%% 
"""
This module is usueful if the data is given as a time series with multiple runs
"""

#Identify the ramps by looking at a min value

min_current = 700 # this is a min value at which all the ramps cross
min_currents_intercepts_1 = np.where(np.abs(current_on_matrix.loc[:,"isrc"]*1000-min_current) < 50) # we fin the points close to 1000 by 50 (amps)

# filter out sections that are too close
min_current_intercepts = []
for i in range(0,len(min_currents_intercepts_1[0])-1):
    i_check = min_currents_intercepts_1[0][i+1]-min_currents_intercepts_1[0][i]
    if i_check > 10: #if the two intercepts are closer than 10 index points, we skip them 
        min_current_intercepts.append(min_currents_intercepts_1[0][i])

min_current_intercepts.append(min_currents_intercepts_1[0][-1])

min_current_start_stop = np.zeros([int(len(min_current_intercepts)-1), 2])

for i in range(0,(int(len(min_current_intercepts)-1))):
    min_current_start_stop[i,0] = min_current_intercepts[i]
    min_current_start_stop[i,1] = min_current_intercepts[i+1]

#save the different runs in a dictionary
all_runs = {}
for i in range(0,len(min_current_start_stop)):
    all_runs[i] = current_on_matrix.iloc[int(min_current_start_stop[i,0]):int(min_current_start_stop[i,1])]

#verify that the data points look good
plt.figure(dpi=1200) 
plt.plot(current_on_matrix.loc[:,"time"],current_on_matrix.loc[:,"isrc"]*1000)
for i in range(0,len(min_current_intercepts)):
    plt.axvline(x = current_on_matrix.loc[:,"time"].iloc[min_current_intercepts[i]], color = 'red', linewidth = 0.5,linestyle='-.')

#%% Graphing multiple shots in one graph - NO Inductance
v_tap = 5
usable_shots = [5,6,10,11,12]
#usable_shots = [0] # for when doing only 1

results_summary = np.zeros([len(usable_shots),3])
j = -1
for i in usable_shots:
    current_on= all_runs[i]["isrc"]*1000
    max_current_idx = current_on[current_on == max(current_on)] #find the index where the current is max
    max_current_idx2 = np.where(current_on == max_current_idx.iloc[0])
    
    time = all_runs[i]["time"].iloc[0:max_current_idx2[0][0]]
    
    
    x_data = current_on.iloc[0:max_current_idx2[0][0]]
    y_data = all_runs[i][channels[v_tap]].iloc[0:max_current_idx2[0][0]]
    
    first_guess = [9000, 1e-6, 22, 1e-9]

    def func_w_R(x, Ic, V_floor, n, R, Vc = (tap_dist[v_tap]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R

    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    j = j+1
    results_summary[j,0] = popt[0]
    results_summary[j,1] = popt[2]
    results_summary[j,2] = popt[3]
    
    fit_x = np.linspace(x_data.iloc[0],x_data.iloc[-1],1000)
    fit_y = func_w_R(fit_x, *popt)   

    #
    fig = plt.figure(1)
    ax = fig.gca()
    plt.rcParams['figure.dpi'] = 1000

    all_files.sort() #I left this here to make sure the legend is in order later on

    fig.suptitle("Current vs Voltage, Ch: " + channels[v_tap],fontsize=15)

    #Formatting the axes
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_axisbelow(True)

    ax.plot(x_data,y_data,color = "black", linewidth = 0.5, label = "raw data")
    ax.plot(fit_x,fit_y, label = "Fit: Ic = " + str(round(popt[0],1)) + "[A], n = " + str(round(popt[2],2)) + "; R =" + str(round(popt[3]*1e9,2))+ " [nOhms]")
    
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 0.5,linestyle='-.')

    ax.plot()
    ax.set_xlabel("Current [A]", fontsize=10)
    ax.set_ylabel("Voltage [V]", fontsize=10)
    #ax.rcParams['figure.dpi'] = 300
    
    
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.00002))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
ax.plot(fit_x,np.full_like(fit_x,tap_dist[v_tap]*1e-6), label = " Vc = " + str(tap_dist[v_tap]) + "uV; Avg Ic = " + str(round(np.mean(results_summary[:,0]),1)) + "; Avg n = " + str(round(np.mean(results_summary[:,1]),1)) , color = "tab:red", linestyle = "--")
ax.legend(fontsize = 5)
#%% Graphing w inductance

v_tap = 5
usable_shots = [1,2,6,7,8]
#usable_shots = [0] # for when doing only 1

for i in usable_shots:
    current_on= all_runs[i]["isrc"]*1000
    max_current_idx = current_on[current_on == max(current_on)] #find the index where the current is max
    max_current_idx2 = np.where(current_on == max_current_idx.iloc[0])
    
    time = all_runs[i]["time"].iloc[0:max_current_idx2[0][0]]
    x_data = current_on.iloc[0:max_current_idx2[0][0]]
    y_data = all_runs[i][channels[v_tap]].iloc[0:max_current_idx2[0][0]]
    def I_ramp_rate(x, m, b):
        return x*m + b

    first_guess_L = [0.1, 1000]

    popt1, pcov1 = curve_fit(I_ramp_rate,time,x_data,p0=first_guess_L)
    
    first_guess = [9000, 1e-6, 20, 1e-9, 1e-3]
    
    def func_w_R_n_l(x, Ic, V_floor, n, R, L, m = popt1[0], b = popt1[1], Vc = (tap_dist[v_tap]*(1e-6))):
        return Vc*((m*x+b)/Ic)**n + V_floor + (m*x+b)*R + m*L

    #Here we fit the fucntion using the least square routine
    popt, pcov = curve_fit(func_w_R_n_l,time,y_data,p0=first_guess)

    #Here we make a linear space and apply it to the fucntion to check our fit against the raw data
    fit_x = np.linspace(time.iloc[0],time.iloc[-1],1000)
    fit_y = func_w_R_n_l(fit_x, *popt)   

    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort() #I left this here to make sure the legend is in order later on

    fig.suptitle("Current vs Voltage, Ch: " + channels[v_tap],fontsize=15)

    #Formatting the axes
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_axisbelow(True)

    ax.plot(time*popt1[0] + popt1[1],y_data,label = "raw data", color = "black", linewidth = 0.5    )
    ax.plot(fit_x*popt1[0] + popt1[1],fit_y, label = "Fit: Ic = " 
            + str(round(popt[0],1)) 
            + " [A], n = " + str(round(popt[2],2)) 
            + ", R = " + str(round(popt[3]*1e9,2))
            + " [nOhms], L = " + str(round(popt[4]*1e2,2)) + " [uH]"
            + ", Ramp rate = " + str(round(popt1[0],0)) + " [A/s]")
    
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 1,linestyle='-.')

    ax.plot()
    ax.set_xlabel("Current [A]", fontsize=10)
    ax.set_ylabel("Voltage [V]", fontsize=10)
    #ax.rcParams['figure.dpi'] = 300
    ax.legend(fontsize = 5)

    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #ax.show()
ax.plot(fit_x*popt1[0] + popt1[1],np.full_like(fit_x,tap_dist[v_tap]*1e-6), label = " Vc = " + str(tap_dist[v_tap]) +
        " [uV], Ramp rate = " + str(round(popt1[0],0)) + " [A/s]", color = "tab:red", linestyle = "--")



#%%inductiance  - trial 1, provided the r
def I_ramp_rate(x, m, b):
    return x*m + b

first_guess_L = [0.1, 100]

popt1, pcov1 = curve_fit(I_ramp_rate,time,x_data,p0=first_guess_L)

fit_x = np.linspace(time.iloc[0],time.iloc[-1],1000)
fit_y = I_ramp_rate(fit_x, *popt1)

plt.plot(time,x_data, label = "raw")
plt.plot(fit_x,fit_y, label = "fit, m = "+ str(round(popt1[0],1)))
plt.legend()
plt.xlabel("time", fontsize = 10)
plt.ylabel("current", fontsize = 10)
plt.rcParams['figure.dpi']= 400
