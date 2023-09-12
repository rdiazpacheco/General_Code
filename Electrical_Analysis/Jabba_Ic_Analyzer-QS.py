# -*- coding: utf-8 -*-
"""
Queen Snake Analysis
Created on Wed Nov 23 17:25:30 2022

@author: rdiazpacheco
"""

#%%

import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from QS_Data_Index import *


#% Data import
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
num_files = len(all_files)
All_files = {}
no_ch = 8
for j in range(0,num_files):  
    one_filename = all_files[j]
    
    All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","ChStat","Value","Load Cell [MPa]")
    #All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","VTapFilter","rVal","Load Cell [MPa]")
tot_ch = no_ch

#%%
#"""
Tap_name = {
    1:"Lead 1",
    2:"Weld 1",
    3:"Weld 2",
    4:"Weld 3",
    5:"Lead 2",
    6:"Weld 1 & 2",
    7:"BTW Weld 2 & 3",
    8:"Overall"  
    }
#"""


Tap_dist = {
    1:36.5,
    2:23,
    3:22.5,
    4:22.5,
    5:37,
    6:44.5,
    7:29.5,
    8:97.5
    }
Steps1 = {
    0:[500,1000,1500,2000,2500, 3000, 3100, 3200, 3300, 3400, 3500],
    }

#%% Weld sample

Tap_name = {
    1:"POS Lead 1",
    2:"POS Lead 2",
    3:"NEG Lead 1",
    4:"NEG Lead 1",
    5:"QSC-PVJ-1",
    6:"QSC-PVJ-2",
    7:"QSC_1",
    8:"QSC_2",
    }


Tap_dist = {
    1:41,
    2:41,
    3:41,
    4:41,
    5:27.5,
    6:31.55,
    7:4,
    8:5
    }


#%% Fitting without steps all in the file


all_files.sort()

#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 7
Inv_taps = 1
I_start = 500
num_files = len(all_files)
R_ind = 2599
#End of noise Current Value
N_ind = 550
Ic_P = np.zeros([len(all_files),2])
Mag_f = 1
Ic_RG = 2560
IvsV = 10
first_guess = [2200, 1e-6, 15, 1e-9]
Ic_P = np.zeros([len(all_files),3])
#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
    fit_y = func_w_R(fit_x, *popt)
    
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
    fig.suptitle("Queen Snake Weld Experiment [pre/post] " + str(ch_no) + " - " + Tap_name[ch_no] + "; " + fname[0:8],fontsize=35)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(0.00003))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)
    
    #Cleaned - inductive voltage
    ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float), 
            signal.decimate(Inv_taps*Mag_f*
                            (All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)
                             -Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV),
             linewidth = 1.25, alpha = 0.85-(j/50), linestyle = "-", color = "black") #label = "Raw - V(ind)",
   
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label =  fname[0:8] + " Ch." + str(ch_no)+ ": " + Tap_name[ch_no] 
            + "; Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) 
            + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", 
            linewidth = 3, linestyle = "-.", alpha = 1-(j/50))


    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 25)#, ncols = 4)
       
    #Average at stepped locations
    Ic_P[j,0] = popt[0] #Ic 
    Ic_P[j,1] = popt[2] #nvalue
    Ic_P[j,2] = popt[3]*(1e9) #r value
    

ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", alpha = 1, color = 'r',
        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV"
        +", IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A] - StDev: " + str(round(np.std(Ic_P[:,0]),1)) 
        +", n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)) 
        +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,2]),1)) + " [ n\u03A9]")
    
   
ax.legend(fontsize = 25) #, ncols = 2)


#%% Fitting for joint resistance: Ch 1

all_files.sort()
#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 4
Inv_taps = 1
#Starting current
I_start = 500
#End of resisitve Current Value
num_files = len(all_files)
IvsV = 10
R_ind = 1500
#End of noise Current Value
N_ind = 550
#Mag_f = 1e-5
#IvsV = 1


first_guess = [100e-9, 1e-5]
#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])
    
    x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
    y_data = signal.decimate(Inv_taps*Mag_f*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
    
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
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_xlabel("Current [A]", fontsize=30)
    ax.set_ylabel("Voltage [V]", fontsize=30)   
    #"""
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.yaxis.set_major_locator(MultipleLocator(0.0001))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    
    #"""
    #Plots 
    
    ax.plot(x_data,y_data, linewidth = 1.25, linestyle = "-", color = "black", alpha = 0.95-(j/100))
    ax.plot(fit_x,fit_y,
            label = Tap_name[ch_no] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
            linewidth = 2.5, linestyle = "--", color = "tab:red")
    
     
    
    ax.legend(fontsize = 30)




#%% Fitting Ic, n,  R stepped
"""
Tap_name = {
    1:"QS-CW PVJ 1",
    2:"QS-CW PVJ 2",
    3:"QS-CW PVJ 3",
    4:"NOT CONNECTED",
    5:"POS Lead A",
    6:"Pos Lead B",
    7:"NEG Lead A",
    8:"NEG Lead B"  
    }
#"""

ramp_rates = {
    0:100
}


setup1 = 0
Steps = Steps1[0]

F_start = 0
ch_no = 1
Inv_taps = 1
I_start = 200
num_files = len(all_files)
R_ind = 3499
#End of noise Current Value
N_ind = 300
Ic_P = np.zeros([len(all_files),3])
Mag_f = 1
Ic_RG = 2560
IvsV = 10



for j in range(0, len(all_files)):    
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,600,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy() 
    #y-data
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))
    
    first_guess = [2300, 1e-6, 9, 1e-9]
    def func_w_R(x, Ic, V_floor, n, R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)   

    plt.rcParams["figure.figsize"] = [25, 15]
    #fig, ax = plt.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    all_files.sort()
    fname1 = (all_files[File_num].partition('\\')[2])  #----------- name here

    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Queen Snake Experiment 5 [pre/post-weld] Ch:" + str(ch_no) + " - " + Tap_name[ch_no] + "; " + fname[0:8],fontsize=35)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)

    #Raw
    ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
            signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
            linewidth = 0.75, linestyle = "-", color = "black")
    
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = fname[0:8] + ": Fit, Ch."  +  str(ch_no) + "; Fit: n=" + str(round(popt[2],1)) + ", Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ",
            linewidth = 2.5, linestyle = "-.", alpha = 0.7)
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    

    #Average at stepped locations
    Ic_P[j,0] = popt[0] #Ic 
    Ic_P[j,1] = popt[2] #nvalue
    Ic_P[j,2] = popt[3]*(1e9) #r value
    

ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = "red", linewidth = 5,
           label = "Fit data, Ch.:" + str(ch_no) 
           +", IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A]"
           +", n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)) 
           +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,2]),1)) + " [ n\u03A9]"
           , zorder = 4)  
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
  
ax.legend(fontsize = 30)


#%% Fitting R stepped
Tap_name = {
    1:"QS-CW PVJ 1",
    2:"QS-CW PVJ 2",
    3:"QS-CW PVJ 3",
    4:"NOT CONNECTED",
    5:"POS Lead A",
    6:"Pos Lead B",
    7:"NEG Lead A",
    8:"NEG Lead B"  
    }
#"""



ch_no = 1
first_guess = [40e-9, 1e-5]
ramp_rates = {
    0:50,
    1:50,
    2:100
}

Steps1 = {
    0:[200,300,400,500,600,700,800,900,1000],
    1:[200,300,400,500,600,700,800,900,1000],
    2:[200,300,400,500,600,700,800,900,1000]
    }
setup1 = 0
Steps = Steps1[1]
#%%
F_start = 0
ch_no =  1
Inv_taps = 1
I_start = 200
num_files = len(all_files)
R_ind = 999
#End of noise Current Value
N_ind = 300
Ic_P = np.zeros([len(all_files),2])
Mag_f = 1
Ic_RG = 2560
IvsV = 10
first_guess = [100e-9,1e-5]

for j in range(0, len(all_files)):    
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,500,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy()           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))    
    
    #fit function
    popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_only_R(fit_x, *popt)   

    plt.rcParams["figure.figsize"] = [25, 15]
    #fig, ax = plt.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    all_files.sort()
    fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Queen Snake Curved Wedge " + str(ch_no) + " - " + Tap_name[ch_no] + "; " + fname[0:8],fontsize=35)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)

    #Raw
    ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
            signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
            linewidth = 0.75, linestyle = "-", color = "black")
    
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = "Fit, R=" + str(round(popt[0]*1e9,0)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[j]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)

    #Average at stepped locations
    Ic_P[j,0] = popt[0]
    

ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = "red", linewidth = 5, label = "Fit data, Ch.:" + str(ch_no) +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,0]*1e9),1)) + " [nOhms]", zorder = 4)  

  
ax.legend(fontsize = 30)