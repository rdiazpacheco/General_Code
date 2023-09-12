# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:09:46 2023
Analyser for IV data from Graphana
Big Tuna Sample 13 20230823

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

all_data1 = {}
for i in range(0,len(all_files)):
    filename = all_files[i]
    all_data1[i] = pd.read_csv(filename, header=0)#, skiprows = range(7,24))


#%%
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
    1:270.2,
    2:50.2,
    3:75.1,
    4:100,
    5:49.6,
    6:74.5,
    7:49.8,
    }
#%% This will show where the cosntant current relates to the current ramp
x_array  = np.arange(current_array.index[0],current_array.index[-1]+1)
plt.plot(x_array, current_array)    
for i in range(0,len(current_steps)):
    plt.axvline(x = current_steps[i,0], color = 'red', linewidth = 0.5,linestyle='-.')
    plt.axvline(x = current_steps[i,1], color = 'blue', linewidth = 0.5,linestyle='-.')
    #plt.axvline(x = test_matrix.loc[current_steps[i,0],"ISRC"]*1000, color = 'red', linewidth = 0.5,linestyle='-.')
    #plt.axvline(x = test_matrix.loc[current_steps[i,1],"ISRC"]*1000, color = 'blue', linewidth = 0.5,linestyle='-.')
#%% Graphing multiple shots in one graph - NO Inductance
colors = {
    0: "red",
    1: "blue"
    }

v_tap = 7

min_current = 0
target_current = 2900


results_summary = np.zeros([len(all_files),3])
for i in range(0,len(all_files)):
    #we have to find the min and max indices for each run
    I_min_idx = min(np.where(all_data1[i].loc[:,"ISRC"]*1000 > min_current)[0])
    I_max_idx = max(np.where(all_data1[i].loc[:,"ISRC"]*1000 == max(all_data1[i].loc[:,"ISRC"]*1000))[0])
    I_target_idx = min(np.where(all_data1[i].loc[:,"ISRC"]*1000 > target_current)[0])
    #noise_avg = np.absolute(np.mean(all_data1[i].loc[all_data1[i].loc[:,"ISRC"] == 0,channels[v_tap]]))
    noise_avg = np.mean(all_data1[i].loc[all_data1[i].loc[:,"ISRC"] == 0,channels[v_tap]])
    current_array = all_data1[i].loc[:,"ISRC"][I_min_idx:I_max_idx+1]*1000
    #first we find the places where the current stays the same (Constant current)
    current_steps1 = []
    for j in range(min(current_array.index),max(current_array.index)-1,1):
        step = current_array.loc[j+1]-current_array.loc[j]
        if step < 10:
            current_steps1.append(j)

    #then we find where the current jumps
    current_steps2 = []
    current_steps2.append(current_steps1[0])
    for j in range(0,len(current_steps1)-1):
        step = current_steps1[j+1]-current_steps1[j]
        if step > 1:
            #print(i) 
            current_steps2.append(current_steps1[j])
            current_steps2.append(current_steps1[j+1])
    current_steps2.append(current_steps1[-1])

    #This will group the start and end of the ramps and also list what the currents are
    current_steps = np.zeros([int(0.5*len(current_steps2)),4])
    for j in range(0,len(current_steps)):
        current_steps[j,0] = int(current_steps2[2*j])
        current_steps[j,1] = int(current_steps2[2*j+1])
        current_steps[j,2] = all_data1[i].loc[current_steps[j,0],"ISRC"]*1000
        current_steps[j,3] = all_data1[i].loc[current_steps[j,1],"ISRC"]*1000
    
    #x data corresponds to the steps - we can also add the min current
    x_data = current_steps[:,2]
    x_data = np.insert(x_data,0,min_current)
    y_data = np.zeros(len(current_steps)-1)
    for k in range(0,len(current_steps)-1):
        y_data[k] = np.abs(np.mean(all_data1[i].loc[int(current_steps[k,0]):int(current_steps[k,1]),channels[v_tap]]))
    y_data = np.insert(y_data,0,noise_avg)
    y_data = np.append(y_data,all_data1[i].loc[int(current_steps[-1,0]),channels[v_tap]])
    
    
    first_guess = [2000, 1e-6, 22, 1e-9]

    def func_w_R(x, Ic, V_floor, n, R, Vc = (tap_dist[v_tap]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R

    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    #j = j+1
    results_summary[i,0] = popt[0]
    results_summary[i,1] = popt[2]
    results_summary[i,2] = popt[3]
    
    fit_x = np.linspace(x_data[0],x_data[-1],1000)
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

    ax.plot(x_data,y_data,color = "black", linewidth = 0.5)
    ax.plot(fit_x,fit_y, label = "Fit: Ic = " + str(round(popt[0],1)) + "[A], n = " + str(round(popt[2],2)) + "; R =" + str(round(popt[3]*1e9,2))+ " [nOhms]", linewidth = 0.5)
    ax.scatter(x_data,y_data,s=80, facecolors='none', edgecolors=colors[i])
    
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 0.5,linestyle='-.')

    #ax.plot()
    ax.set_xlabel("Current [A]", fontsize=10)
    ax.set_ylabel("Voltage [V]", fontsize=10)
    #ax.rcParams['figure.dpi'] = 300
    
    
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
ax.plot(fit_x,np.full_like(fit_x,tap_dist[v_tap]*1e-6), label = " Vc = " + str(tap_dist[v_tap]) + " uV; Avg Ic = " + str(round(np.mean(results_summary[:,0]),1)) + " A; Avg n = " + str(round(np.mean(results_summary[:,1]),1)) , color = "tab:red", linestyle = "--")
#ax.scatter(x_data,y_data)
ax.plot(x_data,y_data,color = "black", linewidth = 0.5, label = "raw data")
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

#%% Find where the current is greater than 0

""" 
Many files in one folder
1. set a minimum current
2. find the max current
3. set a max current to fit
"""
test_matrix = all_data1[0]

#%%
min_current = 100
target_current = 2900


I_min_idx = min(np.where(test_matrix.loc[:,"ISRC"]*1000 > min_current)[0])
I_max_idx = max(np.where(test_matrix.loc[:,"ISRC"]*1000 == max(test_matrix.loc[:,"ISRC"]*1000))[0])
I_target_idx = min(np.where(test_matrix.loc[:,"ISRC"]*1000 > target_current)[0])
noise_avg = np.mean(test_matrix.loc[test_matrix.loc[:,"ISRC"] == 0,"ISRC"])

#%% find the locations where the current is stagnant
current_array = test_matrix.loc[:,"ISRC"][I_min_idx:I_max_idx]*1000

#first we find the places where the current stays the same (Constant current)
current_steps1 = []
for i in range(min(current_array.index),max(current_array.index)-1,1):
    step = current_array.loc[i+1]-current_array.loc[i]
    if step < 10:
        current_steps1.append(i)

#then we find where the current jumps
current_steps2 = []
current_steps2.append(current_steps1[0])
for i in range(0,len(current_steps1)-1):
    step = current_steps1[i+1]-current_steps1[i]
    if step > 1:
        #print(i) 
        current_steps2.append(current_steps1[i])
        current_steps2.append(current_steps1[i+1])
current_steps2.append(current_steps1[-1])

#This will group the start and end of the ramps and also list what the currents are
current_steps = np.zeros([int(0.5*len(current_steps2)),4])
for i in range(0,len(current_steps)):
    current_steps[i,0] = int(current_steps2[2*i])
    current_steps[i,1] = int(current_steps2[2*i+1])
    current_steps[i,2] = test_matrix.loc[current_steps[i,0],"ISRC"]*1000
    current_steps[i,3] = test_matrix.loc[current_steps[i,1],"ISRC"]*1000
