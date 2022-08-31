# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 21:54:04 2022

@author: rdiazpacheco
"""

import numpy as np

def equation_a(a,b,c):
    f = (3*a)/(np.sqrt(b+c))
    return    f

def eggs_per_serving(no_people):
    cakes = int(no_people/2)
    eggs = 3*cakes
    flour = .5*cakes
    return (no_people,eggs,flour)
#%%

answer = equation_a(2, 3, 4)
print(answer)

#%%

no_cakes = eggs_per_serving(100)
print(no_cakes)

#%% Eggs per person (per cake)

data = []

for i in range(1,100):
    data.append(eggs_per_serving(i+1))

from matplotlib import pyplot as plt

xpoints = np.zeros([len(data),1])
ypoints1 = np.zeros([len(data),1])
ypoints2= np.zeros([len(data),1])
for i in range(0,len(data)):
    xpoints[i] = data[i][0]  #people
    ypoints1[i] = data[i][1] #eggs
    ypoints2[i] = data[i][2] #flour

#%%


plt.plot(xpoints, ypoints1, label = "eggs")
plt.plot(xpoints, ypoints2, label = "flour")
plt.legend( fontsize=15)
plt.show()