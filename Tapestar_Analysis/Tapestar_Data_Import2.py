# -*- coding: utf-8 -*-
"""
Created on Tue May  3 15:15:56 2022

@author: rdiazpacheco
"""

"""
Created on Mon Dec  6 10:15:50 2021

Humidity studies, Tapestar file exraction

Summary: 
    Takes a .dat file (extracted from the viewer) with the Ic values with respect to x values.
    These values are then binned to a desired size turning the continues measurement into samples. 
    The samples are then grpahed and an average value is given. 
    A graph of the samples is saved. 

@author: rdiazpacheco
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import pandas as pd

#File Directory


#os.chdir('CFS Dropbox\Ruben Rui Diaz-Pacheco\Testing\CSMC\Analysis\CSMC-Analysis\Tapestar_Data\TS_Files')

#%%
File_name = 'CSMC_0W_001_002.dat';


#os.chdir('TS_Files')
#pick file

TS_data_raw = np.genfromtxt('CSMC_0W_001_002.dat',skip_header = 2);

#%%

#Graph the Ic
xpoints = TS_data_raw[:,0];
Ic = TS_data_raw[:,1];


#Find where the IC is not 0, a LOW-PASS filter is set 
Ic_HTS = np.asarray(np.where(TS_data_raw[:,1] >100));
Ic_spacers = np.asarray(np.where(TS_data_raw[:,1] <75));

#The array of  HTS Ic has gaps corresponding to the spacers, this function finds these gaps
xdist2 = np.full_like(Ic_HTS,0);
for x in range(0,len(Ic_HTS[0])-1):
    xdist2[0,x] = Ic_HTS[0,(x+1)]-Ic_HTS[0,x];

#When the distance between 2 large Ic is detected, a gap has been found, 
#this function is a derivative of the low-pass filter set above
non_0 = np.asarray(np.where(xdist2 > 7))
non0_xdist = np.unique(xdist2[0,non_0[1,:]]);


#The jumps/gaps in HTS Ic are translated back to the Ic map
Ic_jumps = np.full_like(non0_xdist,0);
for x in range(0,len(non0_xdist)):
    temp_a = np.asarray(np.where(xdist2[0,:] == non0_xdist[x]));
    Ic_jumps[x] = temp_a[0,0];
    print(Ic_jumps[x])
    
    
Ic_HTS_jumps = np.arange(2*len(Ic_jumps))
for x in range(0,len(Ic_jumps)):
    Ic_HTS_jumps[2*x] = Ic_HTS[0,Ic_jumps[x]];
    Ic_HTS_jumps[2*x+1] = Ic_HTS[0,Ic_jumps[x]+1];
    
#The location of the IC jumps is saved here    
Ic_jumps2 = xpoints[Ic_HTS_jumps];
#Plotting for confirmation
ypoints3 = np.full_like(TS_data_raw[:,0],0);
ypoints3[Ic_HTS_jumps] = 150;

fig, ax = plt.subplots(figsize=(12, 6))

x = Ic
y = ypoints3


# Set general font size
plt.rcParams['font.size'] = '16'

# Set tick font size
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
	label.set_fontsize(20)
	
ax.plot(y, color='Red', label='Low Ic Points')
ax.plot(x, label='Sine wave')

plt.xlabel('Distance [mm]', fontsize=20)
plt.ylabel('Critical Current [Amps]', fontsize=20)

fig.suptitle("TapeStar data from 4 samples subject Orbital Welds",fontsize=25)

plt.show()



#%%
for x in range(0,len(non0_xdist)-1):
    temp_a = np.asarray(np.where(xdist2[0,:] == non0_xdist[x]));
    if len(temp_a > 1):
        for y in range(0,len(temp_a[0])-1):
            Ic_jumps[x] = temp_a[0,y]
            x = x+1
    else: 
        Ic_jumps[x] = temp_a[0,0];
        print(Ic_jumps[x])

#%%


plt.plot(, , label = "HTS-connector interface")
plt.plot(xpoints, Ic)
plt.title("TapeStar data from 4 samples subject Orbital Welds",fontsize=20)
plt.xlabel("Distance [mm]",fontsize=20)
plt.ylabel("Critical Current [Amps]",fontsize=20)
# Set tick font size
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
	label.set_fontsize(16) 
plt.show()
#%%
#Rebbinning the samples using the locations/indeces of the HTS gaps
sample_length = 50; #50 mm
n_tapes = int(input('How many tapes?'))

if len(Ic_HTS_jumps) > (2*n_tapes):
    Ic_HTS_jumps3 = np.delete(Ic_HTS_jumps,[0,(len(Ic_HTS_jumps)-1)]);
else:
    Ic_HTS_jumps3 = Ic_HTS_jumps;

#The Tape's IC and position are saved in 2 dictionaries below
tapesx = {}
tapesIc = {}
for x in range(0,n_tapes):
    tapesx["tape{0}".format(x+1)] = xpoints[Ic_HTS_jumps3[2*x]:Ic_HTS_jumps3[2*x+1]];
    tapesIc["tape{0}".format(x+1)] = Ic[Ic_HTS_jumps3[2*x]:Ic_HTS_jumps3[2*x+1]];
    
# %%
    
#Each of the tapes will not be cut into 50 mm samples.
#The position array will be first divided into 50 

samples_per_tape = ((tapesx["tape1"].max() - tapesx["tape1"].min())//sample_length);
saples_per_tapeindx = int((len(tapesx["tape1"]))//samples_per_tape);
    
samples = {}
x=0;
#for x in range(0,n_tapes):
for y in range(0,saples_per_tapeindx):
    samples["sample_{0}".format(y+1),1] = tapesx["tape{0}".format(x+1)][(y)*saples_per_tapeindx:(y+1)*saples_per_tapeindx]
    samples["sample_{0}".format(y+1),2] = tapesIc["tape{0}".format(x+1)][(y)*saples_per_tapeindx:(y+1)*saples_per_tapeindx]

    
    
    tapesIc["tape{0}".format(x+1)] = Ic[Ic_HTS_jumps3[2*x]:Ic_HTS_jumps3[2*x+1]];

# %%












