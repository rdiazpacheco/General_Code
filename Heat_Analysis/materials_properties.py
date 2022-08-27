"""
Supercurrent wand materials properties
last modified 12/22/21
author: JLC 
"""
import pandas as pd
import numpy as np
import os
from scipy.interpolate import interp1d


import Constants


def fit_func(T,coefs):
    """ modified from owen's heat3.py 
    returns thermal conductivity value from the General NIST fit function as a function of temperature
    source data from https://trc.nist.gov/cryogenics/materials/materialproperties.htm 
    Args:
        T: [float] temperature [K]
        coefs: [set of floats] coefficients specified by NIST 
    Returns:
        k: [float] thermal conductivity [W/m-K]
    """
    n = np.arange(len(coefs))[:,np.newaxis]
    coefs = np.array(coefs)[:,np.newaxis]
    t2= [T for val in enumerate(coefs)]
    T=np.array(t2)[:,np.newaxis]
    k= np.power(10,np.sum(coefs*np.log10(T)**n,0))[0]

    return k


def thermal_conductivity_of(material, T):
    """ modified from owen's heat3.py
    returns the thermal conductivity of a material at a given temperature based on NIST material curves 
    source data from https://trc.nist.gov/cryogenics/materials/materialproperties.htm 
    Args:
        material: [string] must be a case in thermal_conductivity_of
        T: [float] temperature [K]
    Returns:
        k: [float] thermal conductivity [W/m-K]
    """
    if material==Constants.MATERIAL.CU_R3_50.value: # # RRR = 50: Copper uses a different fit function than the general one 
        a,b,c,d,e,f,g,h,i = coefs = (1.8743,-0.41538,-0.6018,0.13294,0.26426,-0.0219,-0.051276,0.0014871,0.003723)
        k = np.power(10,(a + c*T**0.5 + e*T + g*T**1.5 + i*T**2)/(1 + b*T**0.5 + d*T + f*T**1.5 + h*T**2))
        
    elif material==Constants.MATERIAL.CU_R3_100.value: # RRR = 100:
        a,b,c,d,e,f,g,h,i = coefs = (2.2154,-0.47461,-0.88068,0.13871,0.29505,-0.02043,-0.04831,0.001281,0.003207)
        k = np.power(10,(a + c*T**0.5 + e*T + g*T**1.5 + i*T**2)/(1 + b*T**0.5 + d*T + f*T**1.5 + h*T**2))

    elif material==Constants.MATERIAL.CU_R3_150.value: # RRR= 150:
        a,b,c,d,e,f,g,h,i = coefs = (2.3797,-0.4918,-0.98615,0.13942,0.30475,-0.019713,-0.046897,0.0011969,0.0029988)        
        k = np.power(10,(a + c*T**0.5 + e*T + g*T**1.5 + i*T**2)/(1 + b*T**0.5 + d*T + f*T**1.5 + h*T**2))

    elif material==Constants.MATERIAL.SS_316.value:
        coefs = (-1.4087,1.3982,0.2543,-0.626,0.2334,0.4256,-0.4658,0.165,-0.0199)
        k = fit_func(T, coefs)

    elif material==Constants.MATERIAL.G_10.value:
        coefs = (-4.1236,13.788,-26.068,26.272,-14.663,4.4954,-0.6905,0.0397,0)
        k = fit_func(T, coefs)
        
    elif material=="Solder":  # Assumes a .1mm layer. From Ekin p. 63
        k = np.ones_like(T) * 40 *.001/.01**2
        
    elif material=="Sapphire":
        
        k = fit_func(T, coefs)
    else:
        return np.zeros_like(T)

    return k


def specific_heat_of(material,T):
    """ returns the specific heat of a material at a given temperature based on NIST material curves 
    source data from https://trc.nist.gov/cryogenics/materials/materialproperties.htm 
    Args:
        material: [string] must be a case in thermal_conductivity_of
        T: [float] temperature [K]
    Returns:
        cp: [float] Specific heat [J/g-K]
    """
    if (material==Constants.MATERIAL.CU_R3_50.value or material==Constants.MATERIAL.CU_R3_100.value or material==Constants.MATERIAL.CU_R3_150.value):  
        coefs = (-1.91844,-0.15973,8.61013,-18.996,21.9661,-12.7328,3.54322,-0.3797,0)

    elif material==Constants.MATERIAL.SS_316.value:
        if T<51: # data range of equation applicability per NIST table
            coefs = (12.2486,-80.6422,218.743,-308.854,239.5296,-89.9982,3.15315,8.44996,-1.91368)
        else:
            coefs = (-1879.464,3643.198,76.70125,-6176.028,7437.6247,-4305.7217,1382.4627,-237.22704,17.05262)

    elif material==Constants.MATERIAL.G_10.value:
        coefs = (-2.4083,7.6006,-8.2982,7.3301,-4.2386,1.4294,-0.24396,0.015236,0)

    else:
        return np.zeros_like(T)
    cp = fit_func(T, coefs)

    return cp/Constants.kG_TO_G    


def resistivity_of(material,T):
    """ returns the resistivity of a material at a given temperature  
    Args:
        material: [string] must be a case in Resistivity_of
        T: [float] temperature [K]
    Returns:
        rho: [float] resistivity Ohm*m

    Data from Appendix 6.5a,b in Ekin p.575-576, in Ohm*m*10^-8    
    Wiedemann - Franz - Lorenz : thermal conductivity = LN * T / resistivity
    """
    LN = 2.44*10**-8# V^2 / K^2
    # LN 2.23 # copper
    # rho=LN*T/thermal_conductivity_of(material,T)

    #print("Aliya")

    T_points = [10,20,50,77,100,150,200,250,295]
    if material ==Constants.MATERIAL.CU_R3_100.value: 
        R_points =  [0.015,0.017,0.084,0.21,0.34,0.70,1.07,1.41,1.70]
        func = interp1d(T_points,R_points,kind='linear',bounds_error=None)
        rho = func(T)*Constants.uOHM_CM_TO_OHM_M
        #rho = func(T)*1e-8
    elif material ==Constants.MATERIAL.SS_316.value:
        R_points = [53.9,53.9,54.9,56.8,58.8,63.8,68.9,73.3,77.1]
        func = interp1d(T_points,R_points,kind='linear',bounds_error=None)
        rho = func(T)*Constants.uOHM_CM_TO_OHM_M
    elif material == Constants.MATERIAL.HTS.value: 
        # source: "Optimization Method for Extracting Stabilizer Geometry and Properties of REBCO Tapes" Riva 2021
        if T<80:
            rho = 0
        else:
            rho= 0.3957*Constants.uOHM_CM_TO_OHM_M
    else: 
        R_points = 1e20 * np.ones_like(T_points)
        func = interp1d(T_points,R_points,kind='linear',bounds_error=None)
        rho = func(T)*Constants.uOHM_CM_TO_OHM_M
   

    return rho


def properties_He(t,p):
    """returns the specific heat and density of helium at a given temperature and pressure based on NIST data saved for different
    isobaric cases in the NIST data subfolder. Data source:
    # https://webbook.nist.gov/cgi/fluid.cgi?ID=C7440597&TUnit=K&PUnit=atm&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page
    Args:
        t: [float] temperature [K]
        p: [float] regulator set over-pressure [psi]
    Returns:
        cp: [float] speficic heat [J/(g-K)]
        density: [float] density [g/l]
        viscosity: [float] viscosity [g/(m-s)]
        k: [float] thermal conductivity [W/(m-K)]
    """
    if p>10 or p<5:
        raise ValueError("input must be between 5 and 10 psi")
    elif t<4 or t>300:
        raise ValueError("input must be between 4 and 300 K")
    else:
        df =pd.read_csv(os.getcwd()+"/NIST data/"+Constants.HE_PROPERTIES_PATH_BASE+str(p)+".txt",delimiter="\t") 
        rounded=round(t*2)/2 # round to nearest 0.5 K

        row=df.loc[df["Temperature (K)"]==rounded]
        density=float(row["Density (mol/l)"])*Constants.HE_MOL_TO_G # g/l
        cp=float(row["Cp (J/mol*K)"])/Constants.HE_MOL_TO_G # J/g-K
        viscosity=float(row["Viscosity (uPa*s)"]*Constants.uPA_TO_G_PER_M_S2) # g/(m-s)  
        k=float(row["Therm. Cond. (W/m*K)"])

        return cp,density,viscosity,k