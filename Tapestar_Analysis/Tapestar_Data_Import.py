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
#%%
#File Directory

os.chdir('CFS Dropbox\Ruben Rui Diaz-Pacheco\Testing\CSMC\Analysis\CSMC-Analysis\Tapestar_Data\TS_Files')

#%%
#File_name = 'CSMC_0W_001_002.dat';

File_name = 'Frankenstein_export_001.dat';
#os.chdir('TS_Files')
#pick file

TS_data_raw = np.genfromtxt('Frankenstein_export_001.dat',skip_header = 2);
#TS_data_raw = np.genfromtxt('CSMC_0W_001_002.dat',skip_header = 2);

#Graph the Ic
xpoints = TS_data_raw[:,0];
Ic = TS_data_raw[:,1];x

#%%



#Find where the IC is not 0, a LOW-PASS filter is set 
HTS_on = [i for i in Ic if i >= 80] 
HTS_off = [i for i in Ic if i <= 20]
Ic_spacers2 = np.asarray(np.where(TS_data_raw[:,1] <20));
Ic_spacers = Ic_spacers2[0]

#del Ic_spacers2, Ic_HTS2

#The array of  HTS Ic has gaps corresponding to the spacers, this function finds these gaps
xdist2 = np.full_like(Ic_HTS,0);
for x in range(0,len(Ic_HTS)-1):
    xdist2[x] = Ic_HTS[x+1]-Ic_HTS[x];

xdist2 = list(xdist2)    
#%%
#distancearrayx = np.linspace(1,len(xdist2),len(xdist2))
plt.plot(xpoints[:], Ic)
plt.title("TapeStar data from 4 samples subject Orbital Welds")
plt.xlabel("Distance [mm]")
plt.ylabel("Critical Current [Amps]")
plt.show()

#%%
    
#When the distance between 2 large Ic is detected, a gap has been found, 
#this function is a derivative of the low-pass filter set above
non_0 = [i for i in xdist2 if i >= 3]
def get_unique_numbers(numbers):

    list_of_unique_numbers = []

    unique_numbers = set(numbers)

    for number in unique_numbers:
        list_of_unique_numbers.append(number)

    return list_of_unique_numbers
non_0u = get_unique_numbers(non_0)


non_0u2 = set(non_0u)
xx = [i for i, item in enumerate(xdist2) if item in non_0u2]
xx3 = xpoints[xx]


    

    #%%
Ic_HTS_jumps = np.arange(2*len(Ic_jumps))
for x in range(0,len(Ic_jumps)):
    Ic_HTS_jumps[2*x] = Ic_HTS[0,Ic_jumps[x]];
    Ic_HTS_jumps[2*x+1] = Ic_HTS[0,Ic_jumps[x]+1];
    
#The location of the IC jumps is saved here    
Ic_jumps2 = xpoints[Ic_HTS_jumps];

#Plotting for confirmation
ypoints3 = np.full_like(TS_data_raw[:,0],0);
ypoints3[Ic_HTS_jumps] = 150;
plt.plot(xpoints, ypoints3, label = "HTS-connector interface")
plt.plot(xpoints, Ic)
plt.show()


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
    
#Each of the tapes will now be cut into 50 mm samples.
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












