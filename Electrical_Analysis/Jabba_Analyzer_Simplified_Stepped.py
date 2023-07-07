# -*- coding: utf-8 -*-
"""
Jabba Analyzer simplified
Only for Stepped fits

Created on Fri Jul  7 16:39:09 2023
@author: rdiazpacheco

Variables to edit:


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

#%% Example INNPUTS unique to each experiment
# Example: Weld 20230706

#The tap names are listed as dictionary. Each of these are tied to a channel number
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

#Tap distances are always in Centimeters!
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

#The Steps is a dictionary of the current steps
Steps1 = {
    0:[500,1000,1500,2000,2500, 3000, 3100, 3200, 3300, 3400, 3500],
    }

#The current ramp rates are listed here. We assume the steps are the same. If the steps are not, then the we have to write out each step and then modify the loops :( 
ramp_rates = {
    0:100
}

#If a tap is inverted, applya -1 (This could be made as an IF/THEN in the future)
Inv_taps = 1

#This is the starting current, usually more than 150 Amps in the current set-up. With higher current ramp rates, this number may have to go up to 500. BUT if a initial step at 250 A is used, this is not necessary.
I_start = 250

#I_target is the Max current we wish to analyze. The name can change since it is not tied to the functions. This must be adjusted if the voltage at the max current is far above 4xVc. For lead fits it's fine, but one should use the Resitive fit below.
#Also, you should use Imax to check the MAX CURRENT. If I_target is more than this, the code will default to 500 Amps.
I_target = 3499

#N_ind is the current at which we'd like to conside the noise. IF a "noise" step was used, then than value should be used.  
N_ind = 300

#If there is a modeled IC, this would go here
Ic_model = 2560

#The first guess follows: [Ic, Voltage floor, N-value, Resistance]
first_guess = [2200, 1e-6, 15, 1e-9]

#These are for old files
Mag_f = 1
IvsV = 10

#%% Fitting and Graphing

#The channel you want to look at:
ch_no = 5

#This is the number of points we want to trim from the edges of the current steps. You will need more points if the Current Ramp rate is large. Adjust it so that the scatter points are centered near the bottom of the "voltage wells"
trim_points = 600

#These here are always the same, I just leave them because old data is messed up. 
all_files.sort()
F_start = 0
setup1 = 0


#If the steps are different in the runs, then we'd modify this part and insert it in the loop below. 
Steps = Steps1[0]

#This is an empty table where we will store the Ic, n and R values
Ic_P = np.zeros([len(all_files),3])

#The loop here only goes through all the files, but only on one channel. One could go through all the channels and all the files easily, but I like to make sure each fit to each channel is ok
for j in range(0, len(all_files)):    
    File_num = j+F_start
    
    #This makes a map of current values to indices and also finds the max current in the run
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    
    #This will find the voltage signals within our determined noise range
    I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
    
    #This will find the voltage signals within our determined initial current and target max current
    I_indices_target = range_between_two_Ivalues(All_files[int(File_num)],I_start, I_target)
    
    #This finds the average signal in the noise range 
    Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
    
    #This will find the average value of the voltage at each of the steps while cutting the edges set by the trim_points variable
    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,trim_points,IvsV,tot_ch,Imax)

    #x-data here
    x_data = Steps.copy() 
    #y-data
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))
    
    #This is the function we use. We call it in the loop because we will re-set the Critical Voltage by multiplying the tap distance [in cm] by 1 uV
    def func_w_R(x, Ic, V_floor, n, R, Vc = (Tap_dist[ch_no]*(1e-6))):
        return Vc*(x/Ic)**n + V_floor + x*R
    
    #Here we fit the fucntion using the least square routine
    popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
    
    #Here we make a linear space and apply it to the fucntion to check our fit against the raw data
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_w_R(fit_x, *popt)   

    #A lot of this is just formatting our graphs
    
    #This sets the graph size
    plt.rcParams["figure.figsize"] = [25, 15]
    
    #Setting up the plot
    fig = plt.figure(1)
    ax = fig.gca()
    all_files.sort() #I left this here to make sure the legend is in order later on
    
    #This extracts the name from the folder and uses parts of it to label the plot
    fname1 = (all_files[File_num].partition('\\')[2])  #----------- name here
    fname = fname1[:-4]
    
    #Title of the plot
    fig.suptitle("Queen Snake Experiment 5 [pre/post-weld] Ch:" + str(ch_no) + " - " + Tap_name[ch_no] + "; " + fname[0:8],fontsize=35)
    
    #Formatting the axes
    ax.tick_params(axis='x', labelsize=25)
    ax.tick_params(axis='y', labelsize=25)
    ax.set_axisbelow(True)
    
    #This section formats the major and minor tic marks. I usually disable it while getting the fit right since a wrong fit will crash the computer by placing a ton of tic marks in a plot. 
    #You can disable it/enable it by adding/deleting the hastag symbol ahead of the 3 quotation marks below.
    #"""
    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.yaxis.set_major_locator(MultipleLocator(0.00005))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')
    #"""
    
    #Setting the name to the axes
    ax.set_xlabel("Current [A]", fontsize=25)
    ax.set_ylabel("Voltage [V]", fontsize=25)

    #Here I plot the Raw data with some legend. I ty to plot all of it and you can play with the thickness of the lines to differentiate them.
    ax.plot(All_files[int(File_num)].iloc[int(Current_indices[0]):int(Current_indices[1]),0].astype(float), 
            signal.decimate(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no].astype(float),10),
            linewidth = 0.75, linestyle = "-", color = "black")
    
    #Here we plot the fitted function, one per file. I also pull from POPT the IC, n, and resistive values. 
    ax.plot(fit_x, fit_y,
            label = fname[0:8] + ": Fit, Ch."  +  str(ch_no) + "; Fit: n=" + str(round(popt[2],1)) + ", Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) + " n\u03A9 ",
            linewidth = 2.5, linestyle = "-.", alpha = 0.7)
    ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
    

    #Here I save the fit values to be later averaged
    Ic_P[j,0] = popt[0] #Ic 
    Ic_P[j,1] = popt[2] #nvalue
    Ic_P[j,2] = popt[3]*(1e9) #r value
    
#Now we are outside the loop, we escaped! Yehaa!

#Here I use one scatter since all the files have the same steps. I also use this scatter to add a lot of the average fit values.
ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = "red", linewidth = 5,
           label = "Fit data, Ch.:" + str(ch_no) 
           +", IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A]"
           +", n$_{avg}$= " + str(round(np.mean(Ic_P[:,1]),1)) 
           +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,2]),1)) + " [ n\u03A9]"
           , zorder = 4)  

#Here I graph the Critical voltage value to show that my fit function is legal.
ax.plot(All_files[int(File_num)].iloc[int(I_indices_target[0]):int(I_indices_target[1]),0].astype(float),
        np.full_like(All_files[int(File_num)].iloc[int(I_indices_target[0]):int(I_indices_target[1]),0],Tap_dist[ch_no]*1e-6),
        linewidth = 3, linestyle = ":", alpha = 1, 
        label = "Vc = " + str(round(Tap_dist[ch_no],1)) + "uV")
  
ax.legend(fontsize = 30)


#%% Fitting R stepped

ch_no = 5

#Here I target is something far from the transition
I_target = 2500

#The ramp rates we consider shoudl be below the I_target, you can adjust it by looking at the fit. 
Steps1 = {
    0:[500,1000,1500,2000,2500],
    }

#These here are always the same, I just leave them because old data is messed up. 
all_files.sort()
F_start = 0
setup1 = 0

#If the steps are different in the runs, then we'd modify this part and insert it in the loop below. 
Steps = Steps1[0]

#This is an empty table where we will store the R and Voltage noise values. We only have 2 factors in this fit
Ic_P = np.zeros([len(all_files),2])

#The First guess here only includes Resistance and voltage floor
first_guess = [100e-9,1e-5]

#Everything here works the same as above. Only the function changes
for j in range(0, len(all_files)):    
    File_num = j+F_start
    Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
    I_indices_target = range_between_two_Ivalues(All_files[int(File_num)],I_start, I_target)
    
    Average_signal_at_steps = average_value_at_step(All_files[int(File_num)],Steps,500,IvsV,tot_ch,Imax)
    
    #x-data here
    x_data = Steps.copy()           
    #Y_data here 
    y_data = Average_signal_at_steps.iloc[ch_no-1,:].values
    y_data = list(map(float, y_data))    
    
    #fit function this fit fucntion just finds the slope
    popt, pcov = curve_fit(func_only_R,x_data,y_data,p0=first_guess)
    fit_x = np.linspace(I_start,Imax,len(All_files[int(File_num)].iloc[int(IvsV*Current_indices[0]):int(IvsV*Current_indices[1]),ch_no]))
    fit_y = func_only_R(fit_x, *popt)   

    #Same formatting stuff
    plt.rcParams["figure.figsize"] = [25, 15]
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
            label = "Fit, R=" + str(round(popt[0]*1e9,0)) + ", [nOhms]; Ramp Rate: " + str(ramp_rates[0]) + " [A/s]", linewidth = 2.5, linestyle = "-.", alpha = 0.7)

    #Average at stepped locations
    Ic_P[j,0] = popt[0]
    

ax.scatter(x_data,y_data, s = 200, facecolors = 'none', edgecolors = "red", linewidth = 5, label = "Fit data, Ch.:" + str(ch_no) +", R$_{avg}$= " + str(round(np.mean(Ic_P[:,0]*1e9),1)) + " [nOhms]", zorder = 4)  

  
ax.legend(fontsize = 30)