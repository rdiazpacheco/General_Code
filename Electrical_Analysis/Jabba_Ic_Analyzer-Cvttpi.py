# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 11:16:39 2023

@author: rdiazpacheco

Cavattapi IC Test
"""

#%%

import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from Cvppi_Summary import *

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
    #All_files[str(j+1)].iloc[:,0].interpolate("linear", inplace = True)

tot_ch = no_ch
all_files.sort()
#%% Parameters for a cable
day = 3
Cable = 1
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

first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20, 1e-9]

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
#%% Fitting without steps all in the file
ch_no = 6
I_start = 200
R_ind = 3999

try: 
    ax.cla()
except:
   j = 0
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
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no-1]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
    fit_y = func_w_R(fit_x, *popt)

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
        label = "Vc = " + str(round(Tap_dist[ch_no-1],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)) + " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
    
ax.legend(fontsize = 20)#,ncols = 2)


#%% Fitting for joint resistance: Ch 1

all_files.sort()
F_start = 0
ch_no = 7
I_start  = 450
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


#%% #One Channel at time stepped function
# Ramp information


Steps = ramp_steps[day][2:]
ramp_rates = ramp_rates_steps[day][2:]
I_start = 120
ch_no = 5
N_ind = 500

#try: 
#    ax.cla()
#except:
#    j = 0
for j in range(0, len(all_files)):    
    first_guess = [Ic_search_parameters[Cable][4], 1e-6, 20, 1e-9]
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])

    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,50,IvsV,tot_ch,Imax)
    Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[int(File_num)],IvsV,tot_ch)[ch_no-1]

    #x-data here
    x_data = Steps.copy()
    x_data.insert(0,500)
    #x_data.insert(len(x_data)+1,Imax-1)
           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))
    y_data.insert(0,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],500)),ch_no]-Avg_inductive_V)
    #y_data.insert(1,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],1000)),ch_no])#-Avg_inductive_V)
    #y_data.insert(len(y_data)+1,All_files[int(File_num)].iloc[int(IvsV*I_idx(All_files[int(File_num)],Imax-1)),ch_no])#-Avg_inductive_V)
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no-1]*(1e-6))):
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
    fig.suptitle("Stepped Current Ramp: Ch. " + str(ch_no) + " - " + Tap_name[ch_no-1] + "; " + fname[0:8], fontsize=35)
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
            label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) +
            " [A], R= " + str(round(popt[3]*1e9,1)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[j]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)

    #Display crical voltage and current
    #ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float),
    #        np.full_like(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no]*1e-6),
    #        linewidth = 3, linestyle = ":", color = "red", alpha = 1)
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    #Average at stepped locations
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = popt[2]
ax.plot(fit_x, fit_y,
        label = "Ch.: " + str(ch_no) + "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) +
        " [A], R= " + str(round(popt[3]*1e9,1)) + ", [nOhms]; " + str(ramp_rates[j]) + " [A/s];" " V$_c$ = " + str(Tap_dist[ch_no-1]) + " [cm]", linewidth = 3, linestyle = "-.", alpha = 0.7, color = Tap_color[ch_no])


ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no-1]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1)

ax.scatter(x_data,y_data, s = 700, facecolors = 'none', edgecolors = Tap_color[ch_no], linewidth = 5, label = "Fit data, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)) + " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%", zorder = 4)  

ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
        signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
        linewidth = 0.75, linestyle = "-", color = "black")
  
ax.legend(fontsize = 20)



#"Ch.: " + str(ch_no) + "; V$_c$ = " + str(Tap_dist[ch_no-1]) + " [cm]"


#%% Effect of Ramp rates
ramp_rates = {
    0:50,
    1:50,
    2:50,
    3:100,
    4:100,
    5:100,
    6:30,
    7:30,
    8:200,
    9:200
    }

Tap_dist2 = {
    1:23,
    2:23.5,
    3:15.5,
    4:10.5,
    5:3.5,
    6:1,
    7:24,
    8:23.5
    }

Tap_color = {
    1: "tab:blue",
    2: "tab:red",
    3: "tab:orange",
    4: "tab:orange",
    5: "tab:pink",
    6: "tab:blue",
    7: "tab:pink",
    8: "tab:blue"
    }

all_files.sort()
F_start = 0
ch_no = 6
Inv_taps = 1
I_start = 500
num_files = len(all_files)
R_ind = 3249
N_ind = 550
Ic_P = np.zeros([len(all_files),2])
Mag_f = 1
Ic_RG = 3080
IvsV = 10

fit_values = []

first_guess = [2900, 1e-6, 20, 1e-9]
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
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_values.append(popt)
    
R_summary = pd.DataFrame(fit_values)
R_summary.rename(columns = {0:"Ic", 1:"Voff", 2:"n", 3:"R"},inplace=True)

#%
aa = []
fig = plt.figure(1)
#fig, ax = plt.subplots()
ax = fig.gca()
#ax.cla()
plt.rcParams["figure.figsize"] = [25, 15]
fig.suptitle("Cavatappi II (Ni Strike) Resistance component vs Distance from leads",fontsize=50)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.set_axisbelow(True)

#"""
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator (1))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
#"""
plt.xlim([0,40])
expon = 1e9
ax.set_xlabel("Distance from leads [cm]", fontsize=35)
ax.set_ylabel("Resitivity component [nOhms]", fontsize=35)
Err_std1 = np.std((R_summary.iloc[:,3]*expon), ddof=1) / np.sqrt(len(R_summary.iloc[0,:]))/Tap_dist[ch_no]
#Err_std1 = np.std((R_summary.iloc[:,3]*expon)/Tap_dist[ch_no], ddof=1) / np.sqrt(len(R_summary.iloc[0,:]))/Tap_dist[ch_no]
Err_std2 = np.std(R_summary.iloc[:,2], ddof=1) / np.sqrt(len(R_summary.iloc[0,:]))

for i in range(0,len(R_summary)-1):
    ax.scatter(Tap_dist2[ch_no], R_summary.iloc[i,3]*expon, c = Tap_color[ch_no], s = 300, marker = "^")
    ax.errorbar(Tap_dist2[ch_no],R_summary.iloc[i,3]*expon,Err_std1, color = 'black', capsize = 10)
i = i+1
ax.scatter(Tap_dist2[ch_no], (R_summary.iloc[i,3]*expon), c = Tap_color[ch_no], s = 300, label = "Ch: " + str(ch_no), marker = "^")
#ax1 = ax.twinx()
    
#for i in range(0,len(R_summary)-1):
#    ax1.scatter(ramp_rates[i], R_summary.iloc[i,2], c = 'tab:red', s = 300, marker = "o")
#    ax1.errorbar(ramp_rates[i],R_summary.iloc[i,2],Err_std2, color = 'black', capsize = 10)
#i = i+1
#ax1.scatter(ramp_rates[i], R_summary.iloc[i,2], c = 'tab:red', s = 300, label = "n value", marker = "o")

#ax1.set_ylabel("n value", fontsize=35)

#ax1.tick_params(axis='x', labelsize=30)
#ax1.tick_params(axis='y', labelsize=30)
    
ax.legend(fontsize = 40)#, loc = 'upper center',bbox_to_anchor=(0.5, 1.))
#ax1.legend(fontsize = 40, loc = 'upper center',bbox_to_anchor=(0.5, 0.9))
    
    



