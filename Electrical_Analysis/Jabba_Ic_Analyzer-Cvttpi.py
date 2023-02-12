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
from QS_Data_Index import *

#%% Data import

folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
num_files = len(all_files)
All_files = {}
no_ch = 8
for j in range(0,num_files):  
    one_filename = all_files[j]
    #All_files[str(one_filename[-25:-4])+"{0}".format(j)] = Extract_voltages_one_file(one_filename,8,4,"CmdAmps.Value","ChStat","Value")
    #All_files[j] = Extract_voltages_one_file(one_filename,no_ch,6,"CmdAmps.Value","ChStat","Value","Load Cell [MPa]")
    All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","ChStat","Value","Load Cell [MPa]")
    #All_files[str(j+1)].iloc[:,0].interpolate("linear", inplace = True)

tot_ch = no_ch

#%% Fitting without steps all in the file
Tap_dist = {
    1:22,
    2:15.5,
    3:20,
    4:71,
    5:40,
    6:15.5,
    7:22,
    8:20
    }
#"""

Tap_name = {
    1:"Cppi_2-Lead 1",
    2:"Cppi_2-Core_1",
    3:"Cppi_2-Core_2",
    4:"Cppi_2-Core_3",
    5:"Cppi_2-Core_4",
    6:"Cppi_2-Core_5",
    7:"Cppi_2_Lead 2",
    8:"Cppi_2-Core_6"
    }
#"""


Tap_color = {
    1: "tab:blue",
    2: "tab:red",
    3: "tab:red",
    4: "tab:orange",
    5: "tab:pink",
    6: "tab:blue",
    7: "tab:blue",
    8: "tab:blue"
    }

all_files.sort()

#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 3
Inv_taps = 1

I_start = 300
num_files = len(all_files)
R_ind = 3249
#End of noise Current Value
N_ind = 500
Ic_P = np.zeros([len(all_files),2])
Mag_f = 1
Ic_RG = 3080
IvsV = 10

first_guess = [2900, 1e-6, 20, 1e-9]
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
    
    def func_w_R(x, Ic, V_floor, n,R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #fit function
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
    fit_y = func_w_R(fit_x, *popt)

    first_guess = [2900, 1e-6, 15, 100e-9]
    
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
    fig.suptitle("Cavatappi II (Ni)", fontsize = 40) #" (" + fname1[0:8] + ") - Jacket vs Fuzz", fontsize = 40)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
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
    
    #Cleaned - inductive voltage
    ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float), 
            signal.decimate(Inv_taps*Mag_f*
                            (All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)
                             -Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV),
             linewidth = 1.25, alpha = 0.85-(j/50), linestyle = "-", color = "black") #label = "Raw - V(ind)",
   
    #Fitted Function
    ax.plot(fit_x, fit_y,
            label = "Ch." + str(ch_no)+ ": " + Tap_name[ch_no] + "; Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ", 
            linewidth = int(15-ch_no), linestyle = "-.", color = Tap_color[ch_no], alpha = 1-(j/50))


    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    ax.legend(fontsize = 25)#, ncols = 4)
    Ic_P[j,0] = popt[0]
    Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]
    
ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", color = Tap_color[ch_no], alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A]," +  " Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
    
ax.legend(fontsize = 25) #, ncols = 2)


#%% Fitting for joint resistance: Ch 1

all_files.sort()
#F_start = int(all_files[0].partition('\\')[2][-5])
F_start = 0
ch_no = 7
Inv_taps = 1
#Starting current
I_start = 300
#End of resisitve Current Value
num_files = len(all_files)
IvsV = 1
R_ind = 2000
#End of noise Current Value
N_ind = 500
Mag_f = 1
IvsV = 10


first_guess = [40e-9, 1e-5]
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
    ax.yaxis.set_major_locator(MultipleLocator(0.00001))
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
    
    ax.legend(fontsize = 20)


#%%
# Ramp information
Steps = {
    1:[3500,3700,3900,4000,4750],
    2:[3500,3700,3900,4000,4750],
    3:[4000,4200,4300,4500,4750],
    4:[4000,4200,4300,4500,4750],
    5:[4000,4200,4300,4500,4750],
    6:[4000,4150,4300,4450,4750],
    7:[4000,4150,4300,4450,5000],    
    8:[4000,4150,4300,4450,5000]
    }

Steps2 = {
    #8:[4100,4200,4300,4500,4750],
    9:[3500,3750,4000,4250,4500],
    0:[4100,4200,4300,4400,4500],
    2:[4100,4200,4300,4400,4750],
    3:[4100,4200,4300,4400,4750],
    }

Tap_dist = {
    1:120,
    2:120,
    3:20,
    4:1180,
    5:790,
    6:790,
    7:790,
    8:20,    
    }

#%% All files extracted, in a large dictionary, where each file has its own matrix


#%%
# Ramp information
Steps = {
    1:[3500,3700,3900,4000,4750],
    2:[3500,3700,3900,4000,4750],
    3:[4000,4200,4300,4500,4750],
    4:[4000,4200,4300,4500,4750],
    5:[4000,4200,4300,4500,4750],
    6:[4000,4150,4300,4450,4750],
    7:[4000,4150,4300,4450,5000],    
    8:[4000,4150,4300,4450,5000]
    }

Steps2 = {
    #8:[4100,4200,4300,4500,4750],
    9:[3500,3750,4000,4250,4500],
    0:[4100,4200,4300,4400,4500],
    2:[4100,4200,4300,4400,4750],
    3:[4100,4200,4300,4400,4750],
    }

Tap_dist = {
    1:120,
    2:120,
    3:20,
    4:1180,
    5:790,
    6:790,
    7:790,
    8:20,    
    }

tot_ch = no_ch

#%% One Channel at time stepped function
#File number
File_num = 2
#Starting current
I_start = 500
#End of resisitve Current Value
num_files = len(all_files)

#ch_no = 6
R_ind = 2000
#End of noise Current Value
N_ind = 2000

#Select noice range, give it a intial and final current where the noise should be
Current_indices, Imax = find_start_end_ramp_onefile(All_files[str(File_num)],I_start)
#Select a range where the voltage can be considered as noise
I_indices_Noise = range_between_two_Ivalues(All_files[str(File_num)],I_start, N_ind)
I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
#Find the average noise/signal between the given ranges

ch_no = 2
Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)
# Stepped portion 
Average_signal_at_steps = average_value_at_step(All_files[str(File_num)],Steps[File_num],100)
Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1]

first_guess = [(Tap_dist[ch_no]*(1e-6)), 4300, Avg_inductive_V, 20] # First guess

#x-data here
x_data = Steps[File_num].copy()
I_rep = 4750
x_data.insert(0,500)
x_data.insert(1,N_ind)
#===============> if some points noisy
del x_data[4:6] #===============> if some points noisy

#==============> If the channel is railed
#I_rep = 4750 #==============> If the channel is railed
#x_data = np.array(x_data) #==============> If the channel is railed
#x_data[-1] = I_rep #==============> If the channel is railed

#Y_data here 
y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
y_data = list(map(float, y_data))
y_data.insert(0,All_files[str(File_num)].iloc[int(10*I_idx(All_files[str(File_num)],500)),ch_no]-Avg_inductive_V)
y_data.insert(1,All_files[str(File_num)].iloc[int(10*I_idx(All_files[str(File_num)],N_ind)),ch_no]-Avg_inductive_V)

del y_data[4:6] #===============> if some points noisy

#if the chanel is railed:
#y_data = np.array(y_data)
#I_replacement = All_files[str(File_num)].iloc[10*I_idx(All_files[str(File_num)],I_rep),ch_no] - Avg_inductive_V
#y_data[-1] = I_replacement


#fit function
popt, pcov = curve_fit(func,x_data,y_data,p0=first_guess)
fit_x = np.linspace(I_start,Imax,len(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]))
fit_x2 = np.linspace(0,Imax,len(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]))
fit_y = func(fit_x, *popt)

plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()
all_files.sort()
fname1 = (all_files[File_num-1].partition('\\')[2])  #----------- name here

fname = fname1[:-4]
plt.rcParams["figure.figsize"] = [25, 15]
fig.suptitle("Ch. " + str(ch_no) + " - " + fname,fontsize=25)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_axisbelow(True)
"""
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(0.01))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
#"""
ax.set_xlabel("Current [A]", fontsize=25)
ax.set_ylabel("Voltage [V]", fontsize=25)

ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0], 
        signal.decimate(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no],10),
        label = "Raw", linewidth = 5, linestyle = "-", color = "black")
#Raw - inductive voltage
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0], 
        signal.decimate(All_files[str(File_num)].iloc[int(10*Current_indices[0]):int(10*Current_indices[1]),ch_no]-Avg_at_NoiseRange_per_tap.iloc[ch_no-1],10),
        label = "Raw - V(ind)", linewidth = 2, linestyle = "-", color = "tab:blue")
#Fitted Function
ax.plot(fit_x, fit_y,label = "Fit, n=" + str(round(popt[3],1)) + " Ic=" + str(round(popt[1],1)), linewidth = 7, linestyle = "-.", color = "tab:red", alpha = 0.7)

#Display crical voltage and current
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],
        np.full_like(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],popt[0]),
        linewidth = 3, linestyle = "-.", color = "red", alpha = 0.7, label = "Fit Vc = " + str(round(popt[0]*1e6,1)) + "uV")
ax.plot(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],
        np.full_like(All_files[str(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", color = "red", alpha = 1, label = "Real Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
ax.axvline(x = round(popt[1],1), color = 'red', linewidth = 2,linestyle='-.')
#Average at stepped locations
ax.scatter(x_data,y_data, s = 700, marker = "v", color = "orange", label = "Fit data", zorder = 4)
ax.legend(fontsize = 20)


#%% Fitting for joint resistance: Ch 1
#%% Fitting R stepped
File_num = 3
ch_no = 1
I_start = 200
R_ind = 2000
I_indices_R = range_between_two_Ivalues(All_files[str(File_num)],I_start, R_ind)
#Avg_at_NoiseRange_per_tap = average_in_range(I_indices_Noise,All_files[str(File_num)],0)
#Avg_at_ResistiveRange_per_tap = average_in_range(I_indices_R,All_files[str(File_num)],0)
#Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-offset_voltage_perCh(All_files[str(File_num)])[ch_no-1]

def func_only_R(x,R,OffSet):
    return x*R+OffSet
first_guess = [40e-9, 1e-5]



x_data = (All_files[str(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0])
y_data = signal.decimate(All_files[str(File_num)].iloc[int(10*I_indices_R[0]):int(10*I_indices_R[1]),ch_no],10)

popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)

fit_x = np.linspace(I_start,R_ind,len(x_data))
fit_y = func_only_R(fit_x, *popt)

all_files.sort()
fname1 = (all_files[File_num-1].partition('\\')[2])
fname = fname1[:-4]
#plt.rcParams["figure.figsize"] = [25, 15]
fig, ax = plt.subplots()
fig.suptitle("Ch. " + str(ch_no) + " - " + fname,fontsize=25)

"""
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.yaxis.offsetText.set_fontsize(20)
ax.set_axisbelow(True)
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(0.00001))
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.set_xlabel("Current [A]", fontsize=25)
ax.set_ylabel("Voltage [uV]", fontsize=25)
"""
#Plots 
ax.plot(x_data,y_data, linewidth = 2, linestyle = "-", color = "black")
ax.plot(fit_x,fit_y,
        label = "R = " + str(round(popt[0]*(1e9),3)) + " n\u03A9 " + "(" + str(round((popt[0]*(1e9))/Tap_dist[ch_no],3)) + " n\u03A9 /cm)",
        linewidth = 5, linestyle = "--", color = "tab:red")
ax.legend(fontsize = 40)

#%
fname = fname + " Ch"+str(ch_no)
plt.savefig(fname + ".png", dpi = 300)



