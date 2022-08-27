"""
Supercurrent wand thermal characteristics calculations 
last modified 12/22/21
author: JLC 

# conductive heat leak + ohmic heating - convective cooling + radiative heat leak (radiative should be negligible = A*emissivity*K*T^4, where K = 5.6*10^-8 W/m^2 K^4)

# heat transfer, first law = m*Cp (T2-T1) --> done

# convective cooling = h*A*(T_surface - T_fluid)

# conductive heat leak = A/L * integral (K(T))*dT [W] --> done

# ohmic heating = P=I*V = I^2 * integral R(T)*dT [W] --> need integral R(T) for OFHC Copper at various RRR if possible
"""

import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
import os
import Constants
import materials_properties
import convective_fluid_dynamics

#%%

def thermal_conductivity_integral(material,t_low,t_high):
    """given a material compute integral under thermal conductivity curve for the temperature gradient along length
    Args:
        material: [string] must be a case in thermal_conductivity_of
        t_low: [float] lower of two temperatures [K]
        t_high: [float] higher of two temperatures [K]
    Returns:
        [float] conductivity integral [W/m]
    """
    return integrate.quad(lambda x: materials_properties.thermal_conductivity_of(material,x),t_low,t_high)[0] # first value is integral, second is error


def compute_conduction_heat_leak(area, length, material, t_low, t_high):
    """given a material and temperature gradient along length of material, return conduction. 
    Args:
        area: [float] cross sectional area [m^2]
        length: [float] length [m]
        material: [string] must be a case in thermal_conductivity_of
        t_low: [float] lower of two temperatures [K]
        t_high: [float] higher of two temperatures [K]
    Returns:
        [float] conduction [W]
    """    
    return area/length*materials_properties.thermal_conductivity_integral(material,t_low,t_high)


def ohmic_heating(material,t,I):
    """ Ohmic heating calculation. 
    Args:
        material: [string] must be a case in resistivity_of
        I_: [float] current [A]
    Returns:
        P: [float] ohmic heating in [W]
    # design requirements : max current 800 A. sample can attain 20 K, preferably also 15 K. current leads do not overheat 
    """
    P= I**2*materials_properties.resistivity_of(material,t)
    return P


def convective_cooling(t_surface,t_fluid,radius,length):
    """ Convective heating calculation.  h*A*(T_surface - T_fluid)
    Args:
        t_surface: [float] temperature at surface of material [K]
        t_fluid: [float] temperature of helium [K]
        radius: [float] radius of probe [m]
        length: [float] length of item of interest [m]
    Returns:
        q: [float] convective cooling in [W]
    """
    h=convective_fluid_dynamics.local_heat_transfer_coefficient(t,radius) # [W/(m^2-K)]
    A = 2*pi*radius*length # surface area [m^2]
    q= h*A*(t_surface-t_fluid) # [W]
    
    return q


def heat_transfer_rate(t_out,t_in, flow_rate, p): 
    """returns the heat transfer rate using first law q = m_dot*cp*(Tout-Tin)
    Args:
        t_out: [float] exiting temperature [K]
        t_in: [float] entering temperature [K]
        flow_rate : [float] flow rate in l/min read from flowmeter
        p: [float] regulator set over-pressure [psi]
    Returns:
        q_dot : heat transfer rate [float] [W]
    """
    t_avg=abs(t_out-t_in)/2 # this is a rough approximation, works best above 20K

    cp,_,_,_=materials_properties.properties_He(t_avg,p) # cp [J/g-K]
    print("cp = "+str(cp) + " J/g-K")

    _,density,_,_=materials_properties.properties_He(Constants.ROOM_TEMP_K,p) # density [g/l]
    print("density at circulation pump = "+str(density) + " g/l")

    m_dot=round(flow_rate/Constants.MIN_TO_SEC*density,Constants.ROUND_PREC) # [g/s]
    print("m_dot = "+str(m_dot) + " g/s")

    q_dot=m_dot*cp*(t_out-t_in)
    print("q_dot = "+str(q_dot) + " W")

    return q_dot


def generate_thermal_conductivity_plot():
    t=np.arange(4.0,298.0,1)
    k50=[materials_properties.thermal_conductivity_of(Constants.MATERIAL.CU_R3_50.value,temp) for temp in t]
    k100=[materials_properties.thermal_conductivity_of(Constants.MATERIAL.CU_R3_100.value,temp) for temp in t]
    k150=[materials_properties.thermal_conductivity_of(Constants.MATERIAL.CU_R3_150.value,temp) for temp in t]
    ss=[materials_properties.thermal_conductivity_of(Constants.MATERIAL.SS_316.value,temp) for temp in t]
    g10=[materials_properties.thermal_conductivity_of(Constants.MATERIAL.G_10.value,temp) for temp in t]

    plt.plot(t,k50,label=Constants.MATERIAL.CU_R3_50.value)
    plt.plot(t,k100,label=Constants.MATERIAL.CU_R3_100.value)
    plt.plot(t,k150,label=Constants.MATERIAL.CU_R3_150.value)
    plt.plot(t,ss,label=Constants.MATERIAL.SS_316.value)
    plt.plot(t,g10,label=Constants.MATERIAL.G_10.value)
    plt.legend()
    plt.title("Thermal Conductivity curves by Temperature")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Thermal Conductivity W/m")
    plt.yscale("log", base=10)
    plt.show()


def generate_conduction_heat_leak_integral_table():
    # compute integrals under thermal conductivity curves and export as table, each integral has 4 K as reference starting temperature
    t=np.arange(4.0,301.0,0.1) # temperature range 4K-300.9 K in 0.1 K steps
    t1=t[0] # 4K
    temperatures=t[1:] # all temps except 4.0 K

    heat_leak50 = [materials_properties.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_50.value,t1,t2) for t2 in enumerate(temperatures)]
    heat_leak100 = [materials_properties.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_100.value,t1,t2) for t2 in enumerate(temperatures)]
    heat_leak150 = [materials_properties.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_150.value,t1,t2) for t2 in enumerate(temperatures)]
    heat_leak_ss = [materials_properties.thermal_conductivity_integral(Constants.MATERIAL.SS_316.value,t1,t2) for t2 in enumerate(temperatures)]
    heat_leak_g10 = [materials_properties.thermal_conductivity_integral(Constants.MATERIAL.G_10.value,t1,t2) for t2 in enumerate(temperatures)]

    df = pd.DataFrame(list(zip(t2,heat_leak50,heat_leak100,heat_leak150,heat_leak_ss,heat_leak_g10)),
        columns=['Temp (K)','OFHC Cu RRR=50 (W/m)','OFHC Cu RRR=100 (W/m)','OFHC Cu RRR=150 (W/m)','316 Stainless (kW/m)','G10 (W/m)'])
    df.to_csv(os.getcwd()+'/Thermal_conductivities_kW_m_integral_table.xlsx', index=False) 






# generate_thermal_conductivity_plot()