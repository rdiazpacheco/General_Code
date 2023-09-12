# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:28:51 2023

Strain Calculator as functions

@author: rdiazpacheco
"""

# Calculate the thermal deformation of the sample
def thermal_deformation(alpha,length,T_initial,T_final):
    dL_thermal = alpha * length * (T_final-T_initial)
    return dL_thermal

def mechanical_deformation_e(length,strain_percent):
    dL_mech = length*(strain_percent/100)
    return dL_mech

def deformation_force_dl(area,length,elastic_modulus,mech_deformation):
    F = mech_deformation * (area * elastic_modulus) / length
    return F

def thermal_strain(alpha,T_initial,T_final):
    e_thermal = alpha * (T_final-T_initial)
    return e_thermal

def mechanical_strain_f(force, area, elastic_modulus):
     strain_mech = force / (area * elastic_modulus)   
     return strain_mech     

def required_length_al(length_sample, total_strain_sample, total_strain_aluminum, total_strain_invar):
    L_aluminum = (length_sample * (total_strain_sample - total_strain_invar)) / (total_strain_aluminum - total_strain_invar)
    return L_aluminum

def required_length_al_2(length_sample, length_invar, total_strain_sample, total_strain_aluminum, total_strain_invar, total_strain_tungsten):
    L_aluminum = (length_sample * (total_strain_sample - total_strain_tungsten)+length_invar * (total_strain_invar - total_strain_tungsten)) / (total_strain_aluminum - total_strain_tungsten)
    return L_aluminum

def required_length_al_3(length_sample, length_tungsten, total_strain_sample, total_strain_aluminum, total_strain_invar, total_strain_tungsten):
    L_aluminum = (length_sample * (total_strain_sample - total_strain_invar)+ length_tungsten * (total_strain_tungsten - total_strain_invar)) / (total_strain_aluminum - total_strain_invar)
    return L_aluminum

def force_from_dL_thermal(thermal_dl_diff, length_sample, length_aluminum, length_invar, area_sample, area_aluminum, area_invar, E_sample, E_aluminum, E_invar):
    force = thermal_dl_diff/((length_aluminum/(area_aluminum*E_aluminum))+(length_sample/(area_sample*E_sample))+(length_invar/(area_invar*E_invar)))
    return force


def force_from_dL_thermal2(thermal_dl_diff, length_sample, length_aluminum, length_invar, length_tungsten, area_sample, area_aluminum, area_invar, area_tungsten, E_sample, E_aluminum, E_invar, E_tungsten):
    force = thermal_dl_diff/((length_aluminum/(area_aluminum*E_aluminum))+(length_sample/(area_sample*E_sample))+(length_invar/(area_invar*E_invar))+(length_tungsten/(area_tungsten*E_tungsten)))
    return force


#def sample_mech_strain()
