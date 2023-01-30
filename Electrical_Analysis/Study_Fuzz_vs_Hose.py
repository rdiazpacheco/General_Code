# -*- coding: utf-8 -*-
"""
Fuzz buttons vs hose clapms study
Created on Wed Jan 25 14:30:16 2023

@author: rdiazpacheco
"""
#%% Dependencies
import os
os.chdir("G:\My Drive\Code\General_Code\Electrical_Analysis")
from Jabba_ParserV2 import *
from QS_Data_Index import *

#%% Data import
folder_path = filedialog.askdirectory()
folder_name = folder_path.partition('IcSearch/')[2]
all_files = glob.glob(os.path.join(folder_path, "*.csv"))
all_files.sort()
num_files = len(all_files)
All_files = {}
no_ch = 8

for j in range(0,12):  
    one_filename = all_files[j]
    All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","VTapFilter","rVal","Load Cell [MPa]")


for j in range(12,num_files):  
    one_filename = all_files[j]
    All_files[j] = Extract_voltages_one_file(all_files,j,one_filename,no_ch,6,"CmdAmps.Value","ChStat","Value","Load Cell [MPa]")
    
tot_ch = no_ch
#%% Looking at one core through Days
QS = QS_A1_1
name_index_1 = {}
for i in range(0,len(all_files)):
    name_index_1[i] = [all_files[i].partition('\\')[2][0:8], all_files[i].partition('\\')[2][-16:-4]]
name_index = pd.DataFrame.from_dict(name_index_1, orient= "index")
del(name_index_1)
name_index = name_index.rename(columns = {0:"Date",1:"Filename"})

ch_matches = {}
for i in range(0,len(QS)):
    array_match = np.where(str(QS[i][0]) == name_index["Date"])
    ch_matches[i] = [str(QS[i][0]), array_match[0]]

#%
#QS = QS_A1_1
tot_ch = no_ch
try: 
    ax.cla()
except:
   i = 0
for i in range(0,len(QS)):
#for i in range(2,3):
    ch_no = QS[i][1]
    tap_dist = QS[i][2]
    tap_type = QS[i][3]
    Inv_taps = QS[i][4]
    I_start = 300
    num_files = len(ch_matches[i][1])
    R_ind = 2350
    N_ind = 500
    IvsV = QS[i][6]
    Mag_f = QS[i][5]
    
    Ic_RG = 2560
    Ic_P = np.zeros([num_files,2])
    first_guess = [2200, 1e-6, 20, 50e-9]
    
    for j in range(0,num_files):
        File_num = ch_matches[i][1][j]
        Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
        I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
        I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
        Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])
        
        x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
        y_data = signal.decimate(Inv_taps*Mag_f*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)-Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
        
        def func_w_R(x, Ic, V_floor, n,R, Vc = (tap_dist*(1e-6))):
            return Vc*(x/Ic)**n + V_floor + x*R
        
        #fit function
        popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
        
        fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
        fit_y = func_w_R(fit_x, *popt)
        
        #Display     
        fig = plt.figure(1)
        ax = fig.gca()
        all_files.sort()
        #fname1 = (all_files[File_num-F_start].partition('\\')[2])  #----------- name here
        
        #fname = fname1[:-4]
        plt.rcParams["figure.figsize"] = [25, 15]
        fig.suptitle("Queen Snake A1",fontsize=40) #+ fname
        ax.tick_params(axis='x', labelsize=30)
        ax.tick_params(axis='y', labelsize=30)
        ax.set_axisbelow(True)
        #"""
        ax.xaxis.set_major_locator(MultipleLocator(500))
        ax.yaxis.set_major_locator(MultipleLocator(0.000025))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(which='major', color='#CCCCCC', linestyle='--')
        ax.grid(which='minor', color='#CCCCCC', linestyle=':')
        #"""
        ax.set_xlabel("Current [A]", fontsize=40)
        ax.set_ylabel("Voltage [V]", fontsize=40)
        
        #Cleaned - inductive voltage
        ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float), 
                signal.decimate(Inv_taps*Mag_f*
                                (All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)
                                 -Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV),
                 linewidth = 0.75, alpha = 0.85-(j/50), linestyle = "-", color = "black") #label = "Raw - V(ind)",
       
        #Fitted Function
        ax.plot(fit_x, fit_y,
                label = "Fit, n=" + str(round(popt[2],1)) + " Ic=" + str(round(popt[0],1)) + " [A], R=" + str(round(popt[3]*(1e9),1)) +
                " n\u03A9, " + str(QS[i][0]) +" "+ QS[i][3], 
                linewidth = 1.5, linestyle = "-.", color = "tab:red", alpha = 1-(j/50))


        ax.axvline(x = round(popt[0],1), color = 'red', linewidth = 2,linestyle='-.')
        ax.legend(fontsize = 25)#, ncols = 4)
        Ic_P[j,0] = popt[0]
        Ic_P[j,1] = All_files[int(File_num)].iloc[0,5]
        
    ax.plot(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float),
            np.full_like(All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0],tap_dist*1e-6),
            linewidth = 3, linestyle = ":", color = "red", alpha = 1, 
            label = "Vc = " + str(round(tap_dist,1)) + "uV, IC$_{avg}$= " + str(round(np.mean(Ic_P[:,0]),1)) + " [A], " +  "Degradaton = " + str(round((100*(1-(np.mean(Ic_P[:,0])/Ic_RG))),2)) + "%")
        
    ax.legend(fontsize = 15, ncol = 3)#, ncols = 4)
        
        
#%% Critical Current / n values / r values
Ic_n_R = []
#QS = QS_A1_1
for i in range(0,len(QS)):
#for i in range(0,2):    
    ch_no = QS[i][1]
    tap_dist = QS[i][2]
    tap_type = QS[i][3]
    Inv_taps = QS[i][4]
    I_start = 300
    num_files = len(ch_matches[i][1])
    R_ind = 2350
    N_ind = 500
    IvsV = QS[i][6]
    Mag_f = QS[i][5]
    first_guess = [2200, 1e-6, 11, 50e-9]
    
    for j in range(0,num_files):
        File_num = ch_matches[i][1][j]
        Current_indices, Imax = find_start_end_ramp_onefile(All_files[int(File_num)],I_start)
        I_indices_Noise = range_between_two_Ivalues(All_files[int(File_num)],I_start, N_ind)
        I_indices_R = range_between_two_Ivalues(All_files[int(File_num)],I_start, R_ind)
        Avg_at_NoiseRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_at_ResistiveRange_per_tap = Mag_f*average_in_range(I_indices_Noise,All_files[int(File_num)],tot_ch,0,IvsV)
        Avg_inductive_V = Avg_at_NoiseRange_per_tap.iloc[ch_no-1]-(Inv_taps*Mag_f*offset_voltage_perCh(All_files[int(File_num)],1,tot_ch)[ch_no-1])        
        x_data = All_files[int(File_num)].iloc[int(I_indices_R[0]):int(I_indices_R[1]),0].astype(float)
        y_data = signal.decimate(Inv_taps*Mag_f*(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no].astype(float)
                                                 -Avg_at_NoiseRange_per_tap.iloc[ch_no-1]),IvsV)
        
        def func_w_R(x, Ic, V_floor, n,R, Vc = (tap_dist*(1e-6))):
            return Vc*(x/Ic)**n + V_floor + x*R        
        #fit function
        popt, pcov = curve_fit(func_w_R,x_data,y_data,p0=first_guess)
        #fit_x = np.linspace(I_start,R_ind,len(All_files[int(File_num)].iloc[int(IvsV*I_indices_R[0]):int(IvsV*I_indices_R[1]),ch_no]))    
        #fit_y = func_w_R(fit_x, *popt)
        Ic_n_R.append([QS[i][0],QS[i][3],*popt])
R_summary = pd.DataFrame(Ic_n_R)
R_summary.rename(columns = {0:"Date", 1:"type", 2:"Ic", 3:"Vfloor", 4:"n", 5:"R"},inplace=True)

#%
aa = []
fig = plt.figure(1)
ax = fig.gca()
ax.cla()
plt.rcParams["figure.figsize"] = [25, 15]
#fig.title("Queen Snake A2",fontsize=40) #+ fname
fig.suptitle("Queen Snake A1 - Resistance vs days",fontsize=50)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.set_axisbelow(True)

#"""
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
#"""
ax.set_xlabel("Date", fontsize=35)
ax.set_ylabel("R nOhms]", fontsize=35)
xx = 5
expon = 1e9
for i in range(0,len(R_summary.iloc[:,0])):
    aa.append(datetime.strptime(str(R_summary.iloc[i,0]), '%Y%m%d').strftime('%m/%d/%Y'))
ab = []
ac = []
for i in range(0,len(aa)):
    #plt.scatter(aa[i],R_summary.iloc[i,2],label = R_summary.iloc[i,1])
    if R_summary.iloc[i,1] == "Fuzz":
        ax.scatter(aa[i],R_summary.iloc[i,xx]*expon, c = 'tab:red', s = 300, marker = "v")
        ab.append(R_summary.iloc[i,xx]*expon)
    else:
        ax.scatter(aa[i],R_summary.iloc[i,xx]*expon, c = 'tab:blue', s = 300, marker = "^")
        ac.append(R_summary.iloc[i,xx]*expon)
Err_std1 = np.std(ab, ddof=1) / np.sqrt(np.size(ab))
Err_std2 = np.std(ac, ddof=1) / np.sqrt(np.size(ac))

ax.scatter(aa[0], R_summary.iloc[0,xx]*expon, c = 'tab:blue', s = 300, label = "Jacket", marker = "^")
ax.scatter(aa[5], R_summary.iloc[5,xx]*expon, c = 'tab:red', s = 300, label = "Fuzz Button", marker = "v")
for i in range(0,len(aa)):
    if R_summary.iloc[i,1] == "Fuzz":
        ax.errorbar(aa[i],R_summary.iloc[i,xx]*expon,Err_std1, color = 'black', capsize = 10)
    else:
        ax.errorbar(aa[i],R_summary.iloc[i,xx]*expon,Err_std2, color = 'black', capsize = 10)
ax.axhline(np.mean(ac), lw = 5, linestyle = ":", color = "tab:blue", label = "Jacket R$_{avg}$= " + str(round(np.mean(ac),2)) + ", Std Err = " + str(round(Err_std2,2)))
ax.axhline(np.mean(ab), lw = 5, linestyle = ":", color = "tab:red", label = "Fuzz R$_{avg}$= " + str(round(np.mean(ab),2))+ ", Std Err = " + str(round(Err_std1,2)))
    
ax.legend(fontsize = 40)





