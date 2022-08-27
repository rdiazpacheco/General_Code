# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:02:10 2020

@author: oduke, edited to remove hard paths - JLC
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

home = str(Path.home())

path = home+"/CFS Dropbox/SPARC_Joint/Research_Projects/RPP008-Tape/Testing/SuperCurrent/SuperCurrent Documentation/system logs/2020-11 System.log"
data = pd.read_csv(path,delimiter='\t',parse_dates=['Time'],na_values='--')
data["deltas"] = data["Time"].diff().dt.total_seconds()

#%%

fine_range = (data["Time"] > pd.Timestamp('2020-11-02 05:57:18')) & (data["Time"] < pd.Timestamp('2020-11-06 21:38:35'))
D = data[fine_range]
# t=20
# range2=((data["Sample (K)"]> t-0.1) & (data["Sample (K)"]<t+0.1))
# D=data[range2]
pd.set_option('display.max_columns', None)
# print(D.describe())

# traces = ["Sample (K)","Insert (K)","S1 (K)","S2 (K)","C1 (K)","C2 (K)","C3 (K)","R1 (K)", "R2 (K)","R3 (K)","R4 (K)"]
traces = ["R1 (K)","R3 (K)","R4 (K)"]
#traces = ["deltas"]
for trace in traces:
    plt.plot(D['Time'],D[trace],'.')
plt.legend(traces)
#plt.xlim(pd.to_datetime("3/11/2020 8:00"),pd.to_datetime("3/11/2020 21:00"))
plt.xlabel("Time")
plt.ylabel("Temperature [K]")
plt.title("Mochi Log 11-02-2020 to 11-06-2020")
plt.show()


#%%
# times = np.array([10*60+50,10*60+57,7*60+31,5*60+23,6*60+36])
# counts = np.array([      9,      10,      8,      6,      6])
# print(times/counts)
# print(np.mean(times/counts))