# -*- coding: utf-8 -*-
"""
JABBA FUNCTIONS
Created on Fri Oct 28 10:06:11 2022

@author: rdiazpacheco
"""

#% Functions

def butter_lowpass_filter(cutoff, fs, order,data):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    #y = []
    y = filtfilt(b, a, data)
    return y

def Extract_Voltages(all_files,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue):
    li = []
    num_files = len(all_files)
    for j in range(0,num_files):        
        df = pd.read_csv(all_files[j], header=header_no, skiprows = range(7,24))
        fname1 = (all_files[j].partition('\\')[2])
        fname = fname1[:-4]
        I_val = df.loc[:,I_col_name];
        li.append(I_val)
        li[((1+no_ch)*j)].rename("I_"+fname,inplace=True)
        for i in range(1,(1+no_ch)):
            VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
            li.append(VTap)
            li[(i)+((1+no_ch)*j)].rename("VTap"+str(i)+"_"+fname,inplace=True)            
        frame = pd.concat(li, axis=1, ignore_index= False)        
    return frame

def Extract_voltages_one_file(filename,no_ch,header_no,I_col_name,V_cols_name,rvalORvalue):
    li = []
    df = pd.read_csv(filename, header=header_no, skiprows = range(7,24))
    fname1 = (all_files[j].partition('\\')[2])
    fname = fname1[:-4]
    I_val = df.loc[:,I_col_name];
    li.append(I_val)
    li[0].rename("I_"+fname,inplace=True)
    for i in range(1,(1+no_ch)):
        VTap = df.loc[:,str(V_cols_name) + "[" +str(i) + "]." + str(rvalORvalue)];
        li.append(VTap)
        li[i].rename("V"+str(i),inplace=True)            
    frame = pd.concat(li, axis=1, ignore_index= False)
    return frame

def get_average_of_tap(files_name,data,Current_indices,tap_of_interest):
    li = [] 
    li.append(pd.Series(data.iloc[:,0].iloc[int(Current_indices[0,0]):int(Current_indices[0,1])].values))

    li1 = []
    for j in range(0,len(files_name)):
        li1.append(data.iloc[:,(1+no_ch)*j+(tap_of_interest)].iloc[int(Current_indices[j,0]):int(Current_indices[j,1])])
        li2 = np.mean(li1, axis = 0)
    li.append(pd.Series(li2))
       
    avg_V_tap = pd.concat(li, axis=1, ignore_index= False)    
    
    return avg_V_tap

def noise_range(files_name,data,Istart,Iend):
    I_nA = []
    I_nB = []   

    I_ABs = np.zeros([len(files_name),2])
    for i in range(0,len(files_name)):
        I_nAp = np.where(data.iloc[:,(1+no_ch)*i] > Istart)
        I_nA.append(I_nAp[0])   
        I_ABs[i,0] = min(I_nA[i])
        #if the current upper bound is bigger than the max current, then we set the noise range 500 steps ahead
        if Iend >= max(InVs_1.iloc[:,(1+no_ch)*i]):  
            I_nB.append(min(I_nA[i]) + 500)
            I_ABs[i,1] = min(I_nA[i]) + 500
        else:
            I_nBp = np.where(data.iloc[:,(1+no_ch)*i] > Iend)
            I_nB.append(I_nBp[0])    
            I_ABs[i,1] = min(I_nB[i])          
    return I_ABs

def find_start_end_ramp(data,Istart):
    I_ind = []
    max_Iall = []
    I_str_stp = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nparray = np.where(data.iloc[:,(no_ch+1)*i]>Istart)
        I_ind.append(I_nparray[0])
        I_max = max(data.iloc[:,(no_ch+1)*i])
        max_Ip = np.where(data.iloc[:,(no_ch+1)*i] == max(data.iloc[:,(no_ch+1)*i]))
        max_Iall.append(max_Ip[0])
        I_str_stp[i,0] = min(I_ind[i])
        I_str_stp[i,1] = min(max_Iall[i])        
    return I_str_stp, I_max

def find_noise_per_tap(Noice_indices,data,filt_before):  
    Avg_noise1 = []
    if filt_before == 1:    
        for j in range(0,len(all_files)):
            for i in range(0,no_ch):
                Avg_noise1.append(np.mean(
                    butter_lowpass_filter(cutoff, fs, order,data.iloc[int(Noice_indices[j,0]):int(Noice_indices[j,1]),j*(1+no_ch)+i+1])))   
    else:
        for j in range(0,len(all_files)):
            for i in range(0,no_ch):
                Avg_noise1.append(np.mean(data.iloc[int(Noice_indices[j,0]):int(Noice_indices[j,1]),j*(1+no_ch)+i+1]))
    Avg_noise = pd.Series(Avg_noise1)   
    if len(all_files) > 1:
        avg_noise_p_tap = []           
        for i in range(0,len(all_files)):
            avg_noise_p_tap.append(np.mean(Avg_noise[[i,i+4,i+8,i+12]]))
    else:
        avg_noise_p_tap = []    
        avg_noise_p_tap.append(np.mean(Avg_noise[i]))                    
    return Avg_noise, avg_noise_p_tap

def average_between_current_values(data, Istart, Iend):#, spacing):
    I_end = []
    I_str = []
    I_str_stp = np.zeros([len(all_files),2])
    for i in range(0,len(all_files)):
        I_nparray = np.where(data.iloc[:,(no_ch+1)*i] == Istart)
        I_end.append(I_nparray[0])
        max_Ip = np.where(data.iloc[:,(no_ch+1)*i] == Iend)
        I_str.append(max_Ip[0])
        I_str_stp[i,0] = min(I_end[i])
        I_str_stp[i,1] = min(I_str[i])    
    
    return I_str_stp

def average_at_current(data,I_select,tap_of_choice):#,spacing):

