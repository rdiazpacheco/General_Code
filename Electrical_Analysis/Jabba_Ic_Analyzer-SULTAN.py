# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:55:21 2023

@author: rdiazpacheco
"""

import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from Sultan_summary import *

setup1 = 2

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
tot_ch = no_ch

#%% Fitting without steps all in the file
#"""
all_files.sort()
F_start = 0
ch_no = 1
Inv_taps = 1
I_start = 250
num_files = len(all_files)
R_ind = 3999
N_ind = 500
Ic_P = np.zeros([len(all_files),2])
Mag_f = 1
IvsV = 10
first_guess = [3999, 1e-6, 10, 1e-9]
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
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no][setup1]*(1e-6))):
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
    fig.suptitle(Tap_name[ch_no][setup1] + " " + fname, fontsize = 40) #" (" + fname1[0:8] + ") - Jacket vs Fuzz", fontsize = 40)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    """
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
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
            label = "Ch." + str(ch_no)+ ": " + Tap_name[ch_no][setup1] + "; Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", 
            linewidth = 3, linestyle = "-.", color = Tap_color[ch_no], alpha = 1-(j/50))


    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 25)#, ncols = 4)
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]
    
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no][setup1]*1e-6),
        linewidth = 3, linestyle = ":", color = Tap_color[ch_no], alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no][setup1],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A]," ) #+  " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
    
ax.legend(fontsize = 25) #, ncols = 2)


#%% Fitting for joint resistance: Ch 1

all_files.sort()
#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 8
Inv_taps = 1
#Starting current
I_start = 300
#End of resisitve Current Value
num_files = len(all_files)
IvsV = 10
R_ind = 2000
#End of noise Current Value
N_ind = 400
#Mag_f = 1e-5
#IvsV = 1


first_guess = [40e-9, 1e-5]
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
    fig.suptitle(Tap_name[ch_no][setup1] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=40)
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
            label = Tap_name[ch_no][setup1] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no][setup1],3)) + " n\u03A9 /cm)",
            linewidth = 2.5, linestyle = "--", color = "tab:red")
    
    Ic_P[j,0] = popt[0]*1e9
    Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]


ax.plot(fit_x,fit_y,
        label = "Average R = "+ str(round(np.mean(Ic_P[:,0]),2)) + " n\u03A9 ", linewidth = 2.5, linestyle = "--", color = "tab:red")   
ax.legend(fontsize = 20)
    
    
#%% #One Channel at time stepped function
# Ramp information
ramp_steps = {
    1:[200, 300, 400, 500, 600, 700, 800, 900],
    2:[200, 300, 400, 500, 600, 700, 800, 900],
    3:[200, 300, 400, 500, 600, 700, 800, 900]
    }


Steps = ramp_steps[1]
ramp_rates = 100
I_start = 200
ch_no = 4
Ic_P = np.zeros([len(all_files),3])

#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):    
    first_guess = [2500, 1e-6, 6, 10e-9]
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,10 ,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy()
           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))

    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no][setup1]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)   
    
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = popt[2]
    Ic_P[j,2] = popt[3]

    plt.rcParams["figure.figsize"] = [25, 15]
    #fig, ax = plt.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    all_files.sort()
    fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("SULTAN Foxtrot: Stepped Current Ramp: Ch. " + str(ch_no) + " - " + Tap_name[ch_no][setup1] + "; " + fname[0:8],fontsize=35)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.0001))
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
            label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + 
            " [A], R= " + str(round(popt[3]*1e9,1)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[j+1]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)
    
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')

        
ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no][setup1]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1)

ax.scatter(x_data,y_data, s = 700, facecolors = 'none', edgecolors = Tap_color[ch_no], linewidth = 5, 
           label = "Fit Data, Ch. " + str(ch_no) + ": Ic$_{avg}$ =" + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$ = "+ str(round(np.mean(Ic_P[:,1]),1)) + ", R$_{avg}$ = " + str(round(np.mean(Ic_P[:,2]*1e9),1)) + " [nOhms]; V$_c$ = " + str(Tap_dist[ch_no][setup1]) + " [cm]", zorder = 4)  

  
ax.legend(fontsize = 25)

#%% Resistance in STEPS

ch_no = 8
first_guess = [40e-9, 1e-5]
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

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,1500,IvsV,tot_ch,Imax)
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
    fig.suptitle("SULTAN Foxtrot: Stepped Current Ramp: Ch. " + str(ch_no) + " - " + Tap_name[ch_no][setup1] + "; " + fname[0:8],fontsize=35)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
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
            label = "Fit, R=" + str(round(popt[0]*1e9,1)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[j+1]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)

    #Average at stepped locations
    Ic_P[j,0] = popt[0]
    

ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = Tap_color[ch_no], linewidth = 5, label = "Fit data, Ch.:" + str(ch_no) +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,0]*1e9),1)) + " [nOhms]", zorder = 4)  

  
ax.legend(fontsize = 30)

#%% Resistance in STEPS
setup2 = 1


first_guess = [40e-9, 1e-5]

channels = {
    0:[3,1],
    1:[4,2],
    2:[5,3],
    3:[6,4],
    4:[3,5],
    }
for k in range(0, len(channels)):
    ch_no = channels[k][setup2]

    for j in range(0, len(all_files)):    
        File_num = j+F_start
        Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
        I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
        I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
        Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])
    
        Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,1500,IvsV,tot_ch,Imax)
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
        fig.suptitle("SULTAN Foxtrot: Stepped Current Ramp: Ch. " + "; " + fname[0:8],fontsize=35)
        ax.tick_params(axis='x', labelsize=25)
        ax.tick_params(axis='y', labelsize=25)
        ax.set_axisbelow(True)
        #"""
        ax.xaxis.set_major_locator(MultipleLocator(1000))
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
        ax.plot(fit_x, fit_y,linewidth = 2.5, linestyle = "-.", c= Tap_color[ch_no], alpha = 0.7)#,
                #label = "Fit, R=" + str(round(popt[0]*1e9,1)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[j+1]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)
    
        #Average at stepped locations
        Ic_P[j,0] = popt[0]
        
    
    ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = Tap_color[ch_no], linewidth = 5, label = Tap_name[ch_no][setup1] + "; Ch:" + str(ch_no) +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,0]*1e9),1)) + " [nOhms]", zorder = 4)  

  
ax.legend(fontsize = 20)