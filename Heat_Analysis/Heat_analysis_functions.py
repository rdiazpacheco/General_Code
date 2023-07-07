# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 21:07:38 2023

@author: rdiazpacheco

Thermal Analysis Busbar Fucntions
"""

#%% Dpenedencies




#%% Equations
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

def mm2m(mms):
    meters = mms/1000
    return meters
    
def Resisitivity_simple(material,length,width,thickness):
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



