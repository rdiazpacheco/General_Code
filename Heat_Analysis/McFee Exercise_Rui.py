# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 12:28:36 2022

@author: rdiazpacheco

McFee 
"""
#%%
"part one - constant K and omega"
from numpy.lib.type_check import _asfarray_dispatcher
#import heat_analysis
import pandas as pd
import math as math
#import Constants
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
array#mport materials_properties
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.interpolate import UnivariateSpline



def Cross_section(I,TL,L):
    Area_1 = np.zeros(len(L));
    Cross_1 = np.zeros(len(L));
    for i in list(range(0,len(L))):
        Area_1[i] = L[i]*I*(((2*K_cu*s_cu)*(TH-TL))**-0.5)
        Cross_1[i] = round((2*((Area_1[i]/math.pi)**0.5)),3)
        CS = Cross_1
    return CS


K_cu = 3.86 #W/cm K Thermal conductivity 
s_cu = 58.7e4 #1/cmOhm conducticity, omega
TH = 300
#%%
#s = cm/ohms (resistivity)
#K = thermal conductivity
#TH = high temp of step
#TL = low temp of step


#L/A = 1/I[2Ks(TH-TL)]^1/2
#Qdot = 

#L/A = (1/I^2)[(s1-s2)*Q_dot1 + (s2-s3)*Q_dot2 + ... + (s_n-1-s_n)*Q_dotn-1+s_n*Q_dot_n]
#Q_dot = I[2*(k/s)_avg(TH-TL)]^.5

#use owens code to determine RMS of steady state ohmic heating power during run

#LI/A = [(s1-s2)*(2(k/s)_avg*(TH-TL))^0.5 + ...]

#%% First Analysis

#rho = 1/s
#R = rho(L/A) --> A rho(L/R)
#P_elec = R*I^2--> R = P_elec/I^2


a
L = 50 #50 cm
rho_cu = 1/s_cu #cm/ohms
Width = 8 #cm

I_range = np.linspace(1,16000,1000)
P_range = np.linspace(1,20000,10000)
L_range = np.linspace(10,100,20)
Th_range = np.linspace(1,5,50)

def Area_Power_Current(L,I,P):
    A_1app = np.zeros(len(I))
    for i in list(range(0,len(I))):
        A_1app[i] = rho_cu*(L/( P/(I[i]**2)))
    return A_1app

def Power_by_thickness(L,I,W,Th):
    P_out = np.zeros(len(I))
    for i in list(range(0,len(I))):
        P_out[i] = (rho_cu*L*I[i]**2)/(Th*W)
                
    return P_out

Areas_1 = Area_Power_Current(50,I_range,1)
Power_out  = Power_by_thickness(50,I_range,Width,2)

fig, sfig = plt.subplots(2)
fig.suptitle('Bus bar')

#sfig[0] = plt.figure()
#plt.subplots_adjust(bottom=0.35)
#sfig[0].suptitle('Required Cross-sectional area per current at different lengths and power', fontsize=20)
#sfig[1].suptitle('Power Dissipated for different thicknesses', fontsize=20)

fig = plt.figure()
plt.subplots_adjust(bottom=0.35)
fig.suptitle('Required Cross-sectional area per current at different lengths and power', fontsize=20)


ax = fig.subplots()
ax.set_xlabel('Current')
ax.set_ylabel('Required Cross-sectional area [cm^2]')
#p = ax.plot(Length_many,CSections_1)
p, = ax.plot(I_range,Areas_1,'o') #This one updates the line

ax_slideP = plt.axes([0.2, 0.2, 0.6, 0.03])
ax_slideL = plt.axes([0.2, 0.1, 0.6, 0.03])
#ax_slideTh = plt.axes([0.2, 0.1, 0.6, 0.03])


P_factor = Slider(ax_slideP, 'Power [Watts]',
                  (min(P_range)), max(P_range), valinit=min(P_range), valstep=5)
P_factor.label.set_size(20)

L_factor = Slider(ax_slideL, 'Length [cm]',
                  (min(L_range)), max(L_range), valinit=50, valstep=1)
L_factor.label.set_size(20)

#Th_factor = Slider(ax_slideL, 'Length [cm]',
                  #(min(Th_range)), max(Th_range), valinit=min(Th_range), valstep=1)
#Th_factor.label.set_size(20)

#A_C_P (L,I,P)
def update(val):
    current_P = P_factor.val
    current_L = L_factor.val
    #current_L = L_factor.val
    Areas_1 = Area_Power_Current(L = current_L, I = I_range, P=current_P)
    p.set_ydata(Areas_1)  #This one updates the line
    ax.set_ylim(0,max(Areas_1))
    #ax.plot(I_range,Areas_1)  
    fig.canvas.draw()

def update2(val):
    current_Th = Th_factor.val
    current_L = L_factor.val
    P_out = Power_by_thickness(L = current_L, I = I_range, W = 8, Th = current_Th)
    p.set_ydata(Areas_1)  #This one updates the line
    ax.set_ylim(0,max(P_out))
    #ax.plot(I_range,Areas_1)  
    fig.canvas.draw()

P_factor.on_changed(update)
L_factor.on_changed(update)
Th_factor.on_changed(update2)
plt.show()


#%% Power w constant RHo and changing dimensions
MinTh = 1
MaxTh = 5
MinWdth = 5
MaxWdth = 15
Min_L = 45
Max_L = 100

L_range = np.linspace(Min_L,Max_L,(10*(Max_L-Min_L))+1)
Th_range = np.linspace(MinTh,MaxTh,(10*(MaxTh-MinTh))+1)
Wdth_range = np.linspace(MinWdth,MaxWdth,(10*(MaxWdth-MinWdth))+1)

def Power_by_thickness(L,I,W,Th):
    P_out = np.zeros(len(I))
    for i in list(range(0,len(I))):
        P_out[i] = (rho_cu*L*I[i]**2)/(Th*W)
                
    return P_out

P_out = Power_by_thickness(50,I_range,W = 7,Th = 2)

fig = plt.figure()
plt.subplots_adjust(bottom=0.4)
fig.suptitle('Power Dissipated VS Current with different Lengths, Widths, and Thicknesses', fontsize=20)
ax = fig.subplots()
ax.set_xlabel('Current',fontsize = 20)
ax.set_ylabel('Power Dissipated [Watts]',fontsize = 20)

#p = ax.plot(Length_many,CSections_1)
p, = ax.plot(I_range,P_out) #This one updates the line
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)

ax_slideL = plt.axes([0.2, 0.3, 0.6, 0.03])
ax_slideTh = plt.axes([0.2, 0.2, 0.6, 0.03])
ax_slideWdth = plt.axes([0.2, 0.1, 0.6, 0.03])

L_factor = Slider(ax_slideL, 'Length [cm]',
                  (min(L_range)), max(L_range), valinit=50, valstep=0.1)
L_factor.label.set_size(20)

Th_factor = Slider(ax_slideTh, 'Thickness [cm]',
                  (min(Th_range)), max(Th_range), valinit=min(Th_range), valstep=0.1)
Th_factor.label.set_size(20)

Wdth_factor = Slider(ax_slideWdth, 'Width [cm]',
                  (min(Wdth_range)), max(Wdth_range), valinit=min(Wdth_range), valstep=0.1)
Wdth_factor.label.set_size(20)

def update2(val):
    current_Th = Th_factor.val
    current_Wdth = Wdth_factor.val
    current_L = L_factor.val
    P_out = Power_by_thickness(L = current_L, I = I_range, W = current_Wdth, Th = current_Th)
    p.set_ydata(P_out)  #This one updates the line
    ax.set_ylim(0,max(P_out))
    #ax.plot(I_range,Areas_1)  
    fig.canvas.draw()
    
L_factor.on_changed(update2)
Th_factor.on_changed(update2)
Wdth_factor.on_changed(update2)
plt.show()

#%% Power with dimensions and different Rho

Len_t1 = 100 #cm
n = 100 #(one subdivision at every cm)
LN_level = 30 #(assume LN2 is at 30 cm)
Temp_coeff = 0.00386 #copper is stil linear at Ln2
Tl = 70



#%%

def Make_Conductor(L,n,Tl):
    Hot_End_resistivity = 1.68e-8 #dependend on material, for copper it is 1/1.68x10^-8
    Hot_End_conductivity = 1/Hot_End_resistivity
    Cold_End_conductivity = Hot_End_conductivity/(1+Temp_coeff*(290-Tl))
    Conduct_subd = (Hot_End_conductivity-Cold_End_conductivity)/n    
    Rod_array = np.zeros((2,n))
    Subd_1 = L/n
    for i in list(range(0,n)):
        Rod_array[0,i] = Subd_1*(i+1)
        Rod_array[1,i] = Hot_End_conductivity + Conduct_subd*(i+1)
    return Rod_array

def Power2_along(cond1,I,Th,W):
    P_out = np.zeros(len(cond1))
    for i in list(range(0,len(cond1))):
        P_out[i] = (con1[i]*L***2)/(Th*W)
                
    return P_out

Array1 = Make_Conductor(50,100,70)
Power2_along(Array1[1,:],5000,5,10)




#%%\\
#Now with moving variables

TL_many = np.linspace(250,290,50) #min 4, max 80
I_many = np.linspace(1,16000,50) #min 1, max 800
Length_many = np.linspace(10,100,100)

CSections_1 = Cross_section(I = 1,TL = 4,L = Length_many);

# Plotting
#plt.ion()
fig = plt.figure()
plt.subplots_adjust(bottom=0.35)
fig.suptitle('Required Cross-sectional area per length per current', fontsize=20)

ax = fig.subplots()
ax.set_xlabel('Length of Lead [cm]')
ax.set_ylabel('Required Cross-sectional area')
#p = ax.plot(Length_many,CSections_1)
p, = ax.plot(Length_many,CSections_1,'o') #This one updates the line

ax_slideI = plt.axes([0.2, 0.2, 0.6, 0.03])
ax_slideTL = plt.axes([0.2, 0.1, 0.6, 0.03])
bb1 = plt.axes([0.85, 0.125, 0.1, 0.08])

I_factor = Slider(ax_slideI, 'Current',
                  (min(I_many)), max(I_many), valinit=min(I_many), valstep=5)
I_factor.label.set_size(20)
TL_factor = Slider(ax_slideTL, 'Low Temperature',
                  (min(TL_many)), max(TL_many), valinit=min(TL_many), valstep=1)
TL_factor.label.set_size(20)
#L_factor = Slider(ax_slide, 'Length',0.1, 6, valinit=6, valstep=0.2)

def update(val):
    current_I = I_factor.val
    current_TL = TL_factor.val
    #current_L = L_factor.val
    CSections_1 = Cross_section(I = current_I, TL = current_TL, L=Length_many)
    p.set_ydata(CSections_1)  #This one updates the line
    ax.set_ylim(0,max(CSections_1))
    #ax.plot(Length_many,CSections_1)  
    fig.canvas.draw()
    
def  Capture(event):
    current_I = I_factor.val
    current_TL = TL_factor.val
    CSections_1 = Cross_section(I = current_I, TL = current_TL, L=Length_many)
    ax.plot(Length_many,CSections_1,'+')
    fig.canvas.draw()

   
I_factor.on_changed(update)
TL_factor.on_changed(update)
bcapture = Button(bb1, 'Capture', color='red', hovercolor='green')
bcapture.on_clicked(Capture)
plt.show()

#%% Variable Rho and K


# We first have to determine the Rho and K of Copper
THot = 77 #Kelvin
TLow = 4 # Kelvin
dT = 1 #Segments of analysis - can be higher by choosing a smaller number

#This one gives Electrical conducticity over temp
import scipy.interpolate as spi
T_points = [10,20,50,77,100,150,200,250,295]
R_points =  [0.015,0.017,0.084,0.21,0.34,0.70,1.07,1.41,1.70]
func = spi.interp1d(T_points,R_points,kind='linear',fill_value="extrapolate")

#This one give tehrmal OCnductivity an array
def thermal_conductivity_of2(T):
    a,b,c,d,e,f,g,h,i = coefs = (1.8743,-0.41538,-0.6018,0.13294,0.26426,-0.0219,-0.051276,0.0014871,0.003723)
    k = np.power(10,(a + c*T**0.5 + e*T + g*T**1.5 + i*T**2)/(1 + b*T**0.5 + d*T + f*T**1.5 + h*T**2))
    return k



#def make array of temperatures and coefficients
No_segments = (THot-TLow)//dT
No_segmentsR = (THot-TLow)%dT #have to make sure this remaines an integer because this will become the size of matrix

if No_segmentsR == 0:
    Array_T1 = np.zeros(((No_segments+1),8))
    for i in list(range(len(Array_T1))):
        Array_T1[i,0] = THot-(i*dT)

else:
    Array_T1 = np.zeros(((No_segments+2),8))
    Array_T1[len(Array_T1)-1,0] = TLow
    for i in list(range(len(Array_T1)-1)):
        Array_T1[i,0] = THot-(i*dT)

for i in list(range(len(Array_T1))): #[W/m-K]
    Array_T1[i,1] = thermal_conductivity_of2(Array_T1[i,0])
    
for i in list(range(len(Array_T1))): #Ohm/M
    Array_T1[i,2] = 1/(func(Array_T1[i,0])*Constants.uOHM_CM_TO_OHM_M)

Rod_Array = pd.DataFrame()

Rod_Array['Temp_Index'] = Array_T1[:,0]
Rod_Array['Thermal Conductivity (W/K*m)'] = Array_T1[:,1]
Rod_Array['Electrical Conductivity (1/ohm*m)'] = Array_T1[:,2]

I_t1 = 1000
LIA = 0
gg = np.zeros((len(Rod_Array),4))
kos_sum = Rod_Array.iat[0,1]/Rod_Array.iat[0,2]
kos_avg = kos_sum
Q_dot = ((2*kos_avg*(THot-TLow))**0.5)
gg[0,0] = kos_sum
gg[0,1] = kos_avg
gg[0,2] = I_t1*Q_dot
gg[0,3] = (1/I_t1**2)*(Q_dot*(Rod_Array.iat[1,2]-Rod_Array.iat[0,2]))

for i in range(1,len(Rod_Array),1):
    kos_sum = Rod_Array.iat[i,1]/Rod_Array.iat[i,2] + kos_sum
    kos_avg = kos_sum/(i+1)
    
    Q_dot0I = I_t1*((2*kos_avg*(THot-TLow))**0.5)
    LIA_1 = (1/I_t1**2)*(Q_dot*(Rod_Array.iat[i,2]-Rod_Array.iat[i-1,2]))
    #LIA = (row1.iat[0,2]-row0.iat[0,2])*(2*kos_avg*(TH-TL))**0.5 + LIA 
    
    gg[i,0] = kos_sum
    gg[i,1] = kos_avg
    gg[i,2] = Q_dot0I
    gg[i,3] = LIA_1
    
   
Len_req = 50 #meters
Area_req = np.zeros(len(gg))
for i in range(0,len(gg),1):
    Area_req[i] = (Len_req*I_t1)/gg[i,3]
    
#
#plt.plot(Array_T1[:,0],gg[:,3])
plt.plot(Array_T1[:,0],gg[:,3])
plt.xlabel("Temperature (K)",fontsize = 20)
plt.ylabel("LIA",fontsize = 20)
plt.show()    
    #j=j+1
    #gg[i,0] = kos_avg
#%%
plt.plot(Rod_Array['Temp_Index'],gg[:,2])

plt.xlabel("Temperature (K)",fontsize = 20)
plt.ylabel("Qdot/I",fontsize = 20)
plt.show()

#%% Now use Rod_array to get LIA
T_full_array = list(range(60, 290, 5)) #T_min = 77, T_max = 290, T_step = 1

K = []
S = []
LIA = []
T = []

K, S = generate_arrays(T_full_array, Constants.MATERIAL.CU_R3_100.value)

df_vals = pd.DataFrame()
generate_dataframe(T_full_array, K, S, df_vals)

#LIA, T = LIA_sweep(df_vals, 290, 60, 5)


TH = 285
TL = 60
n = 5
df = df_vals

j = 1
LIA = 0
gg2 = np.zeros((len(df),2))

for i in range(TH, TL+n, -n):
    #gg[j,0] = i
    kos_sum = 0

    for k in range(TH, TH-n*j, -n):

        row = df.loc[df["Temp_Index"]==k]
        kos_sum = row.iat[0,1]/row.iat[0,2] + kos_sum

    kos_avg = (1/j) * kos_sum
    Q_dot2 = I_t1*(2*kos_avg*(TH-TL))**0.5 
    
    #LIoA
    
    gg2[j,0] = i
    #gg2[j,0] = i
    gg2[j,1] = Q_dot20I
    j = j + 1

gg2[0,1] = gg2[1,1]
gg2[len(gg2)-1,1] = gg2[len(gg2)-2,1]
#%%
plt.plot(df['Temp_Index'],gg2[:,1])
plt.xlabel("Temperature (K)",fontsize = 20)
plt.ylabel("Qdot/I",fontsize = 20)

#%%
plt.plot(np.flip(Array_T1[:,0]),gg[:,2])

plt.xlabel("Temperature (K)",fontsize = 20)
plt.ylabel("Qdot/I",fontsize = 20)
plt.show()

    #print("kos avg is " + str(kos_avg))

    #row0 = df.loc[df["Temp_Index"]==i]
    #print("row 0 is " + str(row0))
    #row1 = df.loc[df["Temp_Index"]==i-n]
    #print("row 1 is " + str(row1))

    #LIA = (row1.iat[0,2]-row0.iat[0,2])*(2*kos_avg*(TH-TL))**0.5 + LIA 
    #



#%%def generate_dataframe(T_full_array, K_array, S_array, df):
    
def calc_LIA_Const(TH, TL, n, df):

    '''
    df is a pandas dataframe generated by the function generate dataframe with columns for temperature step, K, S
    TH is integer, high temperature (T0)
    TL is integer, low temperature (Tf)
    n is integer (number of temperature steps) 
    LIA is a float
    '''
    
    j = 1
    LIA = 0

    for i in range(TH, TL+n, -n):

        kos_sum = 0

        for k in range(TH, TH-n*j, -n):

            row = df.loc[df["Temp_Index"]==k]
            kos_sum = row.iat[0,1]/row.iat[0,2] + kos_sum

        kos_avg = (1/j) * kos_sum

        #print("kos avg is " + str(kos_avg))

        row0 = df.loc[df["Temp_Index"]==i]
        #print("row 0 is " + str(row0))
        row1 = df.loc[df["Temp_Index"]==i-n]
        #print("row 1 is " + str(row1))

        LIA = (row1.iat[0,2]-row0.iat[0,2])*(2*kos_avg*(TH-TL))**0.5 + LIA 
        j = j + 1

    #LIA = LIA/100 #convert to cm to check against mcfee
    return LIA #unit: A/m

def LIA_sweep(df, TH, TL, n):

    '''
    df is pandas dataframe, generated by the function generate_dataframe
    TH is integer, high temperature (T0)
    TL is integer, low temperature (Tf)
    n is integer (number of temperature steps) 
    LIA is a float
    '''

    T_full_array = list(range(TH-n, TL-n, -n))

    LIA_array = []

    for i in range(TH-n, TL-n, -n):
        #print(i)
        LIA_temp = calc_LIA_Const(TH, i, n, df)
        LIA_array.append(LIA_temp)

    T_sweep = list(range(TH, TL, -n))

    return LIA_array, T_sweep

#%%

T_full_array = list(range(60, 300, 5)) #T_min = 77, T_max = 290, T_step = 1

K = []
S = []
LIA = []
T = []

K, S = generate_arrays(T_full_array, Constants.MATERIAL.CU_R3_100.value)

df_vals = pd.DataFrame()
generate_dataframe(T_full_array, K, S, df_vals)
LIA, T = LIA_sweep(df_vals, 290, 60, 5)

plt.plot(T, LIA)
plt.xlabel("Temperature (K)")
plt.ylabel("LI/A")
plt.show()

#%%




#%%

#Make an array where in each temperature the K and S are given
#K = THERMAL conductivity 
#S = ELECTRIC conductivity
 
#K of Copper RRR = 50
Cu_example 1 = thermal_conductivity_of("OFHC Copper RRR=50",200)














dsfv








































#%%




