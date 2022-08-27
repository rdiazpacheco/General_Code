# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 07:05:50 2022

@author: rdiazpacheco
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 13:46:12 2022

@author: rdiazpacheco
"""

#pip3 install numpy, pandas, scipy, matplotlib

##Heating per second 
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.interpolate as spi
from matplotlib.widgets import TextBox
#Constants
#Width_1 = #cm
#Length_1 = #cm
#Thickness_1 = #cm
#I_in = 6000 #amps

#Constants
Cu_R = 1.77e-8 #Ohm/m
Cu_den = 997 #kg/m^3
Cu_cap = 389  #J/KgK

# Equations
def in2m(inches):
    meters = 0.0254 * inches
    return meters

def in2mm(inches):
    meters = 25.4 * inches
    return meters

def m2mm(mtrs):
    mmtrs = mtrs*1000
    return mmtrs
    
def R0(length,width,thickness):
    R0 = Cu_R*length/(width*thickness)
    return R0

def R_P(P,wb,lb):   
    P_points = [5,10,20,30,40,50,60] #N/mm^2 1MPa = 1N/mm^2
    R_points = [5000,3100,1500,900,750,700,775] #uOhms/mm^2
    P_R_func = spi.interp1d(P_points,R_points,kind='quadratic',fill_value="extrapolate")
    R_press_A = P_R_func(P)
    R_press = (R_press_A* wb*lb)*(1e-6)
    return R_press

def Ohmic_power(t,w,l,h,current,RP):
    Q_dot = current**2 *((Cu_R*(l+h)/(w*t))+RP)
    return Q_dot


def dTdt(w,t,current):
    dTdt = (current**2 /(w*t)**2)*(Cu_R/(Cu_den*Cu_cap))
    return dTdt

def dT_wdth(width,thickness,current):
    Delta_T = np.zeros((2,len(thickness)))
    for i in range(0,len(thickness),1):
        Delta_T[i,1] = thickness[i]
        Delta_T[i,0] = (current**2 /(width* thickness[i])**2)*(Cu_R/(Cu_den*Cu_cap))
    return Delta_T

def Hz_cooling(w,MaxDT):
    Watts_sqm = (5.92*MaxDT**(1.25))/((1000*w)**0.25)
    return Watts_sqm

def Vt_cooling(h,MaxDT):
    Watts_sqm = (7.66*MaxDT**(1.25))/((h*1000)**0.25)
    return Watts_sqm

def Total_area(t,w,l,h):
    AreaT = (2*w*l)+(2*t*l)+(2*t*w)
    return AreaT

def Tot_AreaH(t,w,l,h):
    Area_tot = (2*l*w)
    return Area_tot

def Tot_AreaV(t,w,l,h):
    Area_tot = 2*(w*h+t*(h+l-t))
    return Area_tot

def Tot_Hz_cooling(t,w,l,h,MaxDT): #w and thickness in m
    Hz_coolP = Hz_cooling(w,MaxDT)*Tot_AreaH(t,w,l,h)
    return Hz_coolP

def Tot_Vt_cooling(t,w,l,h,MaxDT): #w and thickness in m
    Hz_coolP = Vt_cooling(h,MaxDT)*Tot_AreaV(t,w,l,h)
    return Hz_coolP

def Tot_cooling(t,w,l,h,MaxDT):
    Total_cool = Vt_cooling(h,MaxDT)*Tot_AreaV(t,w,l,h) + Hz_cooling(w,MaxDT)*Tot_AreaH(t,w,l,h)
    return Total_cool

def Rate_of_cooling(t,w,l,h,MaxDT,current,RP):
    Heat_rem = Ohmic_power(t,w,l,h,current,RP) - Tot_cooling(t,w,l,h,MaxDT)
    dtdt = Heat_rem/(Cu_cap*Cu_den*w*t*(l+h))
    return dtdt  

def Rate_of_heating(t,w,l,h,current,RP):
    Heat = Ohmic_power(t,w,l,h,current,RP)
    dtdt = Heat/(Cu_cap*Cu_den*w*t*(l+h)) 
    return dtdt        
                    
#def dtdt_I(l,h,MaxDT,current,P):
#    Heating_temp = np.zeros((len(MaxDT_array),5))
#    for i in range(0,len(MaxDT_array),1):
#        Heating_temp[i,0] = RT+MaxDT_array[i]
#        Heating_temp[i,1] = Rate_of_cooling(in2m(0.5),in2m(3),l,,h,MaxDT[i],current,P)
#        Heating_temp[i,2] = Rate_of_cooling(in2m(0.5),in2m(4),l,lb,h,MaxDT[i],current,P)
#        Heating_temp[i,3] = Rate_of_cooling(in2m(0.5),in2m(5),l,lb,h,MaxDT[i],current,P)
#        Heating_temp[i,4] = Rate_of_cooling(in2m(0.5),in2m(6),l,lb,h,MaxDT[i],current,P)
#    return Heating_temp

def TotalT_I(t,w,l,lb,h,RT,Maxtime,current,P):
    timescale = np.linspace(0,Maxtime,Maxtime)
    Total_Temp = np.zeros((len(timescale),3))
    dT = Rate_of_cooling(t,w,l,lb,h,current,P)
    Total_Temp[0,0] = timescale[0]
    Total_Temp[0,1] = dT
    Total_Temp[0,2] = RT
    for i in range(1,len(timescale)):
        Total_Temp[i,0] = i    
        Total_Temp[i,2] = Total_Temp[i-1,2]+Total_Temp[i-1,1]        
        dT = Rate_of_cooling(t,w,l,h,Total_Temp[i,2],current,R_P)        
        Total_Temp[i,1] = dT
    return Total_Temp

def Total_Heating(array,t,w,l,h,RT,ST,I_Max,Rup,Rdwn,Ht,RP,n):
    Total_Temp = array
    Total_Temp[:,0] = np.linspace(0,((I_Max/Rup)+Ht+(I_Max/Rdwn)),n) 
    dt = (((I_Max)/Rup)+Ht+(I_Max/Rdwn))/n
    Total_Temp[:,0] = np.linspace(0,(((I_Max)/Rup)+Ht+(I_Max/Rdwn)),n) 

    i=0
    while Total_Temp[i,0] < (I_Max/Rup):
        Total_Temp[i,1] = Rup*Total_Temp[i,0]
        i+=1
        if Total_Temp[i,0] < (I_Max/Rup):
            continue
        
    while Total_Temp[i,0] < ((I_Max/Rup)+Ht):
        Total_Temp[i,1] = I_Max
        i+=1
        if Total_Temp[i,0] < ((I_Max/Rup)+Ht):
            continue
    k=0
    while Total_Temp[i,0] < ((I_Max/Rup)+Ht+(I_Max/Rdwn)):
        Total_Temp[i,1] = I_Max-Rdwn*Total_Temp[k,0]
        i+=1    
        k+=1
        if Total_Temp[i,0] < ((I_Max/Rup)+Ht+(I_Max/Rdwn)):
            continue
    
    if ST < RT:
        Total_Temp[0,2] = 0 #Rate of heating    
        Total_Temp[0,3] = ST #Temperature
                        
        i=1
        while Total_Temp[i-1,3] < RT:            dT1 = Rate_of_heating(t,w,l,h,Total_Temp[i-1,1],RP)
            Total_Temp[i,2] = dT1
            Total_Temp[i,3] = Total_Temp[i-1,3]+Total_Temp[i,2]*dt                
            #print(Total_Temp[i,3])
            if Total_Temp[i-1,3] > RT:
                continue
            elif i == n-1:
                continue
            else:
                i+=1
                
        while Total_Temp[i,0] < max(Total_Temp[:,0]):
            dT1 = Rate_of_cooling(t,w,l,h,Total_Temp[i-1,3]-RT,Total_Temp[i,1],RP)                    
            Total_Temp[i,2] = dT1
            Total_Temp[i,3] = Total_Temp[i-1,3]+Total_Temp[i,2]*dt
            if Total_Temp[i,0] > max(Total_Temp[:,0]):
                break
            elif i == n-1:
                break
            else:
                i+=1
                
    else:
        i=1
        while Total_Temp[i,0] < max(Total_Temp[:,0]):
            dT1 = Rate_of_cooling(t,w,l,h,Total_Temp[i-1,3]-RT,Total_Temp[i,1],RP)                    
            Total_Temp[i,2] = dT1
            Total_Temp[i,3] = Total_Temp[i-1,3]+Total_Temp[i,2]*dt
            if Total_Temp[i,0] > max(Total_Temp[:,0]):
                break
            elif i == n-1:
                break
            else:
                i+=1
    return Total_Temp

#%% #Initial constants
#Geometry swivel
#w = in2m(4) #[m]
#t = in2m(0.5) #[m]
#l = 0.55 #[m]
#h = 0.15 #[m]
#--------------------------------------------------------------------------
#Current & Ramping
I_Max = 25000 #[amps]
Rup = 15 #[A/s] ---> change to RI_up
Rdwn = 100 #[A/s] ---> change to RI_DWN
Ht = 20 #[secs]
#--------------------------------------------------------------------------
#Temperatures
RT = 300 #[K] Room temperature
#ST = 200 #[K] Starting temperature
#ST2 = 299 #K

#--------------------------------------------------------------------------
#Pressure
P = 20 #MPa
#RP_1 = R_P(P,wb,lb)
#--------------------------------------------------------------------------
#Detail of analysis (t/n)
n = 100
ST = 280 #[]


ST = 299
w = in2m(8)
t = in2m(2)
l = in2m(24)
h = 0.000001 #[m]
wb = in2m(6)
lb = in2m(2.5) #[m]
RP_1 = R_P(P,wb,lb)

Rup = 150 #[A/s] ---> change to RI_up
Rdwn = 100 #[A/s] ---> change to RI_DWN
Ht = 3600 #[secs]
earray = np.zeros((n,4))
Total_Temp2 = earray
Total_Temp2 = Total_Heating(earray,t,w,l,h,RT,ST,I_Max,Rup,Rdwn,Ht,RP_1,n)
Total_Temp1 = Total_Temp2

#%%

fig, ax = plt.subplots(2, figsize=(25, 15))

#fig.suptitle('Jumper Cables Heat Analysis',fontsize=25)
fig.suptitle('Power Supply Flags Heat Analysis',fontsize=25)

ax[0].set_title("Total Temperature v Time",fontsize=20)
ax[1].set_title("Current Ramp",fontsize=20)

ax[1].set_xlabel("Time [s]", fontsize=15)
ax[0].set_ylabel("Temperature [K]", fontsize=15)
ax[1].set_ylabel("Current [A]", fontsize=15)

plt.subplots_adjust(bottom=0.2,right = 0.75)

earray = np.zeros((n,4))
RoomTemparray = np.zeros(len(Total_Temp2))
RoomTemparray[:] = 300

TTemp1, = ax[0].plot(Total_Temp2[:-1,0],Total_Temp1[:-1,3],label="Main Lead", linewidth = 3)
#TTemp2, = ax[0].plot(Total_Temp2[:-1,0],Total_Temp3[:-1,3],label="Busbar Swivel Receiver", linewidth = 3)
#TTemp3, = ax[0].plot(Total_Temp2[:-1,0],Total_Temp4[:-1,3],label="Busbar Vertical", linewidth = 4)
#TTemp4, = ax[0].plot(Total_Temp2[:-1,0],Total_Temp5[:-1,3],label="Power Supply Flag", linewidth = 3)
RoomTemp = ax[0].plot(Total_Temp2[:-1,0],RoomTemparray[:-1], label="Room Temperature 300K", color='red', linewidth=4, linestyle=':')

I_in, = ax[1].plot(Total_Temp2[:,0],Total_Temp2[:,1],label="Current [A]",linewidth = 4)
leg = ax[0].legend(fontsize=20)
leg = ax[1].legend(fontsize=20)

#plt.show()

#TTemp, = plt.plot(Total_Temp[:,0],Total_Temp[:,3],label="With Circular section")
#I_in, = plt.plot(Total_Temp[:,0],Total_Temp[:,1],label="Current Ramp")    

#All Graphing needs bellow
#Geometry
#t_box = plt.axes([0.55, 0.8, 0.1, 0.05])
#w_box = plt.axes([0.7, 0.8, 0.1, 0.05])
#l_box = plt.axes([0.85, 0.8, 0.1, 0.05])
#lb_box = plt.axes([0.55, 0.72, 0.1, 0.05])
#h_box = plt.axes([0.7, 0.72, 0.1, 0.05])

#Temperatures
RT_box = plt.axes([0.8, 0.8, 0.1, 0.05])
#ST_box = plt.axes([0.7, 0.6, 0.1, 0.05])

#Ramp Profile
I_Max_box = plt.axes([0.8, 0.6, 0.075, 0.05])
Ht_box = plt.axes([0.9, 0.6, 0.075, 0.05])
Rup_box = plt.axes([0.8, 0.5, 0.075, 0.05])
Rdwn_box = plt.axes([0.9, 0.5, 0.075, 0.05])



#Pressure & Bolts
P_box = plt.axes([0.8, 0.4, 0.1, 0.05])
#B_box = plt.axes([0.1, 0.05, 0.8, 0.075])
#Ab_box = plt.axes([0.1, 0.05, 0.8, 0.075])

#Detail
#n_box = plt.axes([0.55, 0.22, 0.1, 0.05])

#t_initial_text = t
#w_initial_text = w
#l_initial_text = l
#lb_initial_text = lb
#h_initial_text = h
RT_initial_text = RT
#ST_initial_text = ST
Rup_initial_text = Rup
Rdwn_initial_text = Rdwn
Ht_initial_text = Ht
P_initial_text = P
#n_initial_text = n
I_Max_initial_text = I_Max

#text_box1 = TextBox(t_box, 'Thickness [m]', initial=t_initial_text)#, hovercolor='0.975')
#text_box2 = TextBox(w_box, 'Width [m]', initial=w_initial_text)#, hovercolor='0.975')
#text_box3 = TextBox(l_box, 'Length [m]', initial=l_initial_text)#, hovercolor='0.975')
#text_box4 = TextBox(lb_box, 'Overlap [m]', initial=lb_initial_text)#, hovercolor='0.975')
#text_box5 = TextBox(h_box, 'Height [m]', initial=h_initial_text)#, hovercolor='0.975')

text_box6 = TextBox(RT_box, 'Room T [K]', initial=RT_initial_text)#, hovercolor='0.975')
#text_box7 = TextBox(ST_box, 'Start T [K]', initial=ST_initial_text)#,hovercolor='0.975')

text_box8 = TextBox(I_Max_box, 'Max Current [Amps]', initial=I_Max_initial_text)#, hovercolor='6000')
text_box9 = TextBox(Rup_box, 'Ramp Up [A/s]', initial=Rup_initial_text)#, hovercolor='0.975')
text_box10 = TextBox(Rdwn_box, 'Ramp Down [A/s]', initial=Rdwn_initial_text)#, hovercolor='0.975')
text_box11 = TextBox(Ht_box, 'Hold [secs]', initial=Ht_initial_text)#, hovercolor='0.975')


text_box12 = TextBox(P_box, 'Pressure [MPa]', initial=P_initial_text)#, hovercolor='0.975')
#text_box = TextBox(B_box, 'Thickness', initial=initial_text)
#text_box = TextBox(Ab_box, 'Thickness', initial=initial_text)

#text_box13 = TextBox(n_box, 'Detail [sec/n]', initial=n_initial_text)#, hovercolor='0.975')

plt.show()
