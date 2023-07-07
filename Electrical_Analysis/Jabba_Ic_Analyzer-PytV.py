# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 15:57:24 2023

@author: rdiazpacheco
Python 5 Joggle test Data Analysis 
"""
import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from PythonV_Summary import *

#% Data import

folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
all_files.sort()
num_files = len(all_files)
All_files = {}
no_ch = 8
for j in range(0,num_files):  
    one_filename = all_files[j]
    All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","ChStat","Value","CmdVolts.Value") 

    #All_files[str(j+1)].iloc[:,0].interpolate("linear", inplace = True)

tot_ch = no_ch

#%% Parameters for a cable
day = 0
Cable = 0
Tap_dist = Tap_distance[day][2:]
Tap_name = Tap_names[day][1:]
ramp_rates = ramp_rates_ls[day][2:]

I_start = Ic_search_parameters[Cable][1]
R_ind = Ic_search_parameters[Cable][3]
Ic_RG = Ic_search_parameters[Cable][4]
N_ind = Ic_search_parameters[Cable][2]
Ic_P = np.zeros([len(all_files),2])
num_files = len(all_files)
F_start = 0
Inv_taps = 1
Mag_f = 1
IvsV = 10

first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20, 0]
#first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20]

Tap_color = {
    1: "tab:blue",
    2: "tab:red",
    3: "tab:red",
    4: "tab:orange",
    5: "tab:pink",
    6: "tab:orange",
    7: "tab:blue",
    8: "tab:blue"
    }

ccolor = {
    1:"blue",
    2:"green",
    3:"red",
    4:"blue",
    5:"red",
    6:"purple"    
    }

#%% #One Channel at time stepped function w R component
# Ramp information

all_files.sort()
Steps = ramp_steps[day][2:]
ramp_rates = ramp_rates_steps[day][2:]
I_start = 120
ch_no = 3
R_ind = 3599

stylel = {
    1:"-.",
    2:"-.",
    3:"--",
    4:"-.",
    5:":",
    6:"-.",
    7:"-.",
    8:"-."
    }
first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20, 15e-6]
#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):    
    #
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,100,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy()

           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))
    
    def func(x, Ic, V_floor, n,  Vc = (Tap_dist[ch_no-1]*(1e-6)),):
        return Vc*(x/Ic)**n + V_floor
    
    def func_w_R(x, Ic, V_floor, n, R, Vc = (Tap_dist[ch_no-1]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    #fit_x2 = np.linspace(0,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)   

    plt.rcParams["figure.figsize"] = [25, 15]
    #fig, ax = plt.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    all_files.sort()
    fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Stepped Current Ramp: Ch. " + str(ch_no) + " - " + fname[0:8],fontsize=25)
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.00015))
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
    Max_p1 = round(All_files[int(File_num)].iloc[0,9],1)
    if Max_p1 < 0:
        Max_p = 0
    else:
        Max_p = Max_p1
    
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = "Ch: " + str(ch_no) + " - " + Tap_name[ch_no-1] +", n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) +
            " [A], R = " + str(round(popt[3]*1e9,1)) + ", [nOhms] Max P = " + str(Max_p) + " MPa", linewidth = 2.5, linestyle = stylel[ch_no], alpha = 1, color = ccolor[j+1]) #

    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = popt[2]
    #Data points
    ax.scatter(x_data,y_data, s = 700, facecolors = 'none', edgecolors = ccolor[j+1], linewidth = 5) #label = "Fit data, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),
    ax.axvline(x = round(popt[0],1), color = "tab:red", linewidth = 2,linestyle='-.')


ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no-1]*1e-6),
        linewidth = 3, linestyle = ":", color = Tap_color[ch_no], alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no-1],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)))

ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
        signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
        linewidth = 0.75, linestyle = "-", color = "black")
  
ax.legend(fontsize = 25)

#%% #One Channel at time stepped function w/o r component
# Ramp information

all_files.sort()
Steps = ramp_steps[day][2:]
ramp_rates = ramp_rates_steps[day][2:]
I_start = 120
ch_no = 2
R_ind = 3500

stylel = {
    1:"-.",
    2:"-.",
    3:"--",
    4:"-.",
    5:":",
    6:"-.",
    7:"-.",
    8:"-."
    }
first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20]
#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):    
    #first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20, 0]
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,100,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy()
    #x_data.insert(0,550)
    #x_data.insert(len(x_data)+1,Imax-1)
           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))
   # y_data.insert(0,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],500)),ch_no]-Avg_inductive_V)
    #y_data.insert(1,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],1000)),ch_no])#-Avg_inductive_V)
    #y_data.insert(len(y_data)+1,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],Imax-1)),ch_no])#-Avg_inductive_V)
    
    def func(x, Ic, V_floor, n,  Vc = (Tap_dist[ch_no-1]*(1e-6)),):
        return Vc*(x/Ic)**n + V_floor
    
    def func_w_R(x, Ic, V_floor, n, R, Vc = (Tap_dist[ch_no-1]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    #fit function
    popt, pcov = curve_fit(func,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    #fit_x2 = np.linspace(0,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func(fit_x, *popt)   

    plt.rcParams["figure.figsize"] = [25, 15]
    #fig, ax = plt.subplots()
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    all_files.sort()
    fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Stepped Current Ramp: Ch. " + str(ch_no) + " - " + fname[0:8],fontsize=25)
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
    Max_p1 = round(All_files[int(File_num)].iloc[0,9],1)
    if Max_p1 < 0:
        Max_p = 0
    else:
        Max_p = Max_p1
    
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = "Ch: " + str(ch_no) + " - " + Tap_name[ch_no-1] +", n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) +
            " [A], Max P = " + str(Max_p) + " MPa", linewidth = 2.5, linestyle = stylel[ch_no], alpha = 1, color = ccolor[j+1]) #R= " + str(round(popt[3]*1e9,1)) + ", [nOhms]

    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = popt[2]
    #Data points
    ax.scatter(x_data,y_data, s = 700, facecolors = 'none', edgecolors = ccolor[j+1], linewidth = 5) #label = "Fit data, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),
    ax.axvline(x = round(popt[0],1), color = "tab:red", linewidth = 2,linestyle='-.')


ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no-1]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1)
    
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no-1]*1e-6),
        linewidth = 3, linestyle = ":", color = Tap_color[ch_no], alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no-1],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)))# + " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")

#ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
#        signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
#        linewidth = 0.75, linestyle = "-", color = "black")
  
ax.legend(fontsize = 25)



#%% Effect of pressure

all_files.sort()
Steps = ramp_steps[day][2:]
ramp_rates = ramp_rates_steps[day][2:]
I_start = 120
ch_no = 3
R_ind = 3500
N_ind = 500

ccolor = {
    0:"blue",
    1:"green",
    2:"tab:orange",
    3:"tab:green",
    4:"tab:red",
    5:"tab:blue",
    6:"tab:purple"
    }

stylel = {
    1:"-.",
    2:"-.",
    3:"--",
    4:"-.",
    5:":",
    6:"-.",
    7:"-.",
    8:"-."
    }
first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20]

channels = {
    0:4,
    1:5,
    2:2,
    3:6,
    4:3
    }
fit_values = []

popt2 = []

for k in range(0, len(channels)):
    ch_no = channels[k]
    R_summary = []
    fit_values = []
    popt2 = []
    
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
        
        Max_p1 = round(All_files[int(File_num)].iloc[0,9],1)
        if Max_p1 < 0:
            Max_p = 0
        else:
            Max_p = Max_p1
        
        def func(x, Ic, V_floor, n,  Vc = (Tap_dist[ch_no-1]*(1e-6)),):
            return Vc*(x/Ic)**n + V_floor
        
        #fit function
        popt, pcov = curve_fit(func,x_data,y_data,p0=first_guess)
        popt2 = np.append(popt,Max_p)
        
        fit_values.append(popt2)
        
    R_summary = pd.DataFrame(fit_values)
    R_summary.rename(columns = {0:"Ic", 1:"Voff", 2:"n", 3:"MaxP"},inplace=True)
    aa = []
    fig = plt.figure(1)
    #fig, ax = plt.subplots()
    ax = fig.gca()
    #ax.cla()
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle("Python V N Value vs Pressure",fontsize=50)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_axisbelow(True)

    #"""
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator (1))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    #plt.ylim([2250,3500])
    item_list = 2
    expon = 1e0
    ax.set_xlabel("Applied Pressure [MPa]", fontsize=35)
    ax.set_ylabel("Critical currnet [A]", fontsize=35)


    for i in range(0,len(R_summary)):
        #ch_no2 = i%4
        ax.scatter(R_summary.iloc[i,3], (R_summary.iloc[i,item_list]*expon), c= ccolor[ch_no], s = 300, marker = "^")
        #ax.errorbar(R_summary.iloc[i,3],R_summary.iloc[i,item_list]*expon,Err_std1, color = 'black', capsize = 10)
    #i = i+1
    ax.scatter(R_summary.iloc[i,3], (R_summary.iloc[i,item_list]*expon), c = ccolor[ch_no], s = 300, label = "Ch: " + str(ch_no), marker = "^")

ax.legend(fontsize = 40,loc = "upper right")

#%% Fitting without steps all in the file
ch_no = 2
I_start = 200
R_ind = 3599
first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20]
#try: 
#    ax.cla()
#except:
#   j = 0
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
    
    def func(x, Ic, V_floor, n, Vc = (Tap_dist[ch_no-1]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor
    
    #fit function
    popt, pcov = curve_fit(func,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
    fit_y = func(fit_x, *popt)

    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
    
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    #fig.suptitle(Tap_name[ch_no] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=25)
    fig.suptitle(Tap_distance[day][1] + ", " + str(Tap_distance[day][0]), fontsize = 40) #" (" + fname1[0:8] + ") - Jacket vs Fuzz", fontsize = 40)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_axisbelow(True)
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
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
            label = "Ch." + str(ch_no)+ ": " + Tap_name[ch_no-1] + "; Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 " + "Ramp Rate: " + str(ramp_rates[j]), 
            linewidth = 2, linestyle = "-.", color = Tap_color[ch_no], alpha = 1-(j/50))


    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 25)#, ncols = 4)
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = popt[2]
    
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no-1]*1e-6),
        linewidth = 3, linestyle = ":", color = Tap_color[ch_no], alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no-1],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)))# + " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
    
ax.legend(fontsize = 30)#,ncols = 2)


#%% Fitting for joint resistance: Ch 1

all_files.sort()
F_start = 0
ch_no = 7
I_start  = 200
R_ind = 2000

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
    
    fig = plt.figure(2)
    
    ax = fig.gca()
    all_files.sort()
    fname1 = (all_files[File_num-F_start].partition('\\')[2])
    fname = fname1[:-4]
    fname = fname1[:-4]
    plt.rcParams["figure.figsize"] = [25, 15]
    fig.suptitle(Tap_name[ch_no-1] + ": Ch. " + str(ch_no) + " - " + fname,fontsize=40)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    ax.set_xlabel("Current [A]", fontsize=30)
    ax.set_ylabel("Voltage [V]", fontsize=30)   
    #"""
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(MultipleLocator(250))
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    
    #"""
    #Plots 
    ax.plot(x_data,y_data, linewidth = 1.25, linestyle = "-", color = "black", alpha = 0.95-(j/100))
    ax.plot(fit_x,fit_y,
            label = Tap_name[ch_no-1] + ": R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
            linewidth = 2.5, linestyle = "--", color = "tab:red")
    
    ax.legend(fontsize = 20)