# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 22:18:00 2022

@author: rdiazpacheco
"""
import numpy as np

def find_initial_speed(v_f,dist,acc):
    v_i = np.sqrt(v_f**2 - 2*acc*dist)
    return v_i

def m_s_from_k_hr(kmhr):
    ms =(kmhr*1000/3600)
    ms_str = str(ms) +' m/s'
    return ms_str

def ms_from_mh(mh):
    ms =(mh*1609/3600)
    ms_str = str(ms)+'m/s'
    return ms_str

def distance_w_SAT(speed,acce,time):
    xdist = (0.5*acce**2*time)+speed*time
    return xdist


def distance_w_IFT(init,final,time):
    xdist = ((init+final)/2)*time
    return xdist
#%%
print(str(find_initial_speed(40, 10, 9))+"m/s")