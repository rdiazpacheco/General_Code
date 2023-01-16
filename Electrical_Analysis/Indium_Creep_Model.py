# -*- coding: utf-8 -*-
"""
Modeling Indium Creep and Potential Energy 
Created on Fri Dec 23 15:25:46 2022

@author: rdiazpacheco
"""

#Dependencies
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy.interpolate as spi
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from numpy import diff
#%%
num_cs = 4
thickness = np.linspace(0,2,1000)
tvf = pd.DataFrame(thickness)
tvf2 = pd.DataFrame(thickness)
tvf.loc[:,0].rename("thickness", inplace = True)
tvf2.loc[:,0].rename("thickness", inplace = True)

aa = 1

for i in range(0,num_cs):
    #force = ((1-(i/7))/((i+10)*tvf.loc[:,0]**2 + (0.005*i)))+(0.01-(i*0.01))
    force1 = ((aa-0.2*i)/(2*tvf.loc[:,0])) - (tvf.loc[:,0])
    force2 = ((aa-0.1*i)**2)/(2*(tvf.loc[:,0])**0.5)
    
    tvf[str("C_" +str(i))] = force1
    tvf2[str("C_" +str(i))] = force2

f_val = 1.65
intersections = np.zeros([tvf.shape[1]-1,2])
for i in range(1,tvf.shape[1]):    
    bigger_thanINDX = []
    bigger_than = []
    bthanArray = np.where(tvf.iloc[:,i] > f_val)    
    bigger_thanINDX.append(bthanArray[0])
    intersections[i-1,0] = int(max(bigger_thanINDX[0]))
    intersections[i-1,1] = tvf.iloc[int(intersections[i-1,0]),i]
    
#% PLotting
Temp = {
        1: "300 K - Sphere",
        2: "355 K - Sphere",
        3: "401.15 K (Epoxy curing temp.)",
        4: "430 K"
            }
Temp2 = {
        1: "300 K - Wire",
        2: "325 K - Wire",
        3: "401.15 K (Epoxy curing temp.)",
        4: "430 K"
            }

fig = plt.figure(1)
ax = fig.gca()
ax.cla()

fname1 = "Conceptual Function of Indium Thickness vs Clamping Force"

plt.rcParams["figure.figsize"] = [25, 15]
fig.suptitle(fname1,fontsize=45)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.set_axisbelow(True)
ax.set_xlim([0,1])
ax.set_ylim([0,4])
#"""
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
#"""
ax.set_xlabel("Thickness (\u03C6) [mm]", fontsize=40)
ax.set_ylabel("Radius [mm]  (Clamping Force \u221D $r^2$)", fontsize=40)

x = np.linspace(0, 1, 10)
number = tvf.shape[1]-2
cmap = plt.get_cmap('inferno')
colors = [cmap(i) for i in np.linspace(0, 1, number)]
#Fitted Function
for i, color in enumerate(colors, start=1):
    ax.plot(tvf.iloc[:,0], tvf.iloc[:,i], color = color, label = "T = " + Temp[i], linewidth = 5, linestyle = "-")
    #ax.plot(tvf.iloc[:,0], tvf2.iloc[:,i], color = color, label = "T = " + Temp2[i], linewidth = 5, linestyle = "-")
    
    ax.axvline(x = tvf.iloc[int(intersections[i-1,0]),0], color = 'tab:blue', linewidth = 3,linestyle='-.')
ax.plot(tvf.iloc[:,0], np.full_like(tvf.iloc[:,0],0.01), linewidth = 5, linestyle = ":", color = "gray", alpha = 1, label = "T = 430 K (In melting point)")
ax.plot(tvf.iloc[:,0], np.full_like(tvf.iloc[:,0],f_val), linewidth = 5, linestyle = ":", color = "red", alpha = 1, label = "~30MPa over 200mm")    
ax.legend(fontsize = 40)

    
#%% Energy of creep
fig = plt.figure(1)
ax = fig.gca()
ax.cla()

fname1 = "Conceptual Function of Potential Energy in the Clamp System"

plt.rcParams["figure.figsize"] = [25, 15]
fig.suptitle(fname1,fontsize=45)
ax.tick_params(axis='x', labelsize=30)
ax.tick_params(axis='y', labelsize=30)
ax.set_axisbelow(True)
ax.set_xlim([0,1])
ax.set_ylim([0,15])
#"""
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(2.5))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.set_yticklabels([])
#"""
ax.set_xlabel("Thickness (\u03C6) [mm]", fontsize=40)
ax.set_ylabel("Energy [kJ]", fontsize=40)

x = np.linspace(0,1,500)
x2 = np.linspace(0,1,499)
for j in range(0,3):
    eps_0 = 1-(j*0.35)    
    y1 = (eps_0**5)*(1-(x/eps_0))*((eps_0-((x**2)/eps_0))/(eps_0*x**2))    
    ax.plot(x,y1, linewidth = 5, label = "$\u03C6_0$ = " + str(round(eps_0,3)) + "mm")
y5 = (eps_0**5)*0.5*(1-(x/eps_0))*((eps_0-((x**2)/eps_0))/(eps_0*x**2))
ax.plot(x,y5, linewidth = 5, label = "$\u03C6_0$ = " + str(round(eps_0,3)) + "mm, $G(400K) < G(300K)$" )
   
eps_0 = 0.3
y1 = (eps_0**5)*(1-(x/eps_0))*((eps_0-((x**2)/eps_0))/(eps_0*x**2)) 
y_int2 = np.where(x < 0.2)
y_int1 = max(y_int2[0])
ax.scatter(0.2,y1[y_int1],linewidth = 25, color = "tab:blue", label = "Compressing from 300\u03BCm to 200\u03BCm")

eps_0 = 1    
y1 = (eps_0**5)*(1-(x/eps_0))*((eps_0-((x**2)/eps_0))/(eps_0*x**2)) 
y_int2 = np.where(x < 0.3)
y_int3 = np.where(x < 0.2)
y_int1 = max(y_int2[0])
y_int2 = max(y_int3[0])
ax.scatter(0.3,y1[y_int1],linewidth = 25, color = "tab:red", label = "Compressing from 1mm to 300\u03BCm")
ax.scatter(0.2,y1[y_int2],linewidth = 25, color = "black", label = "Compressing from 1mm to 200\u03BCm w/o heating")
ax.scatter(0.183,y1[y_int1],linewidth = 25, color = "tab:orange", label = "Thickness if $E_{bvlle}$ = $E_{deformation}$")
ax.scatter(0.156,y1[y_int1],linewidth = 25, color = "tab:purple", label = "Thickness if $E_{bvlle}$ = $E_{deformation}$ &  $G(400K) < G(300K)$")

ax.axvline(x = (0.3), color = 'tab:red', linewidth = 3,linestyle='-.')    
ax.axvline(x = (0.2), color = 'tab:blue', linewidth = 3,linestyle='-.')
ax.axvline(x = (0.183), color = 'tab:orange', linewidth = 3, linestyle='-.')
ax.axvline(x = (0.156), color = 'tab:purple', linewidth = 3, linestyle='-.')
    
ax.legend(fontsize = 20, loc = 1)



    