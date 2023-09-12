# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 00:56:19 2023

@author: rdiazpacheco

Strain calculator

"""
####################################
# Thermal Stress Test V2
# Solve for L_aluminum and L_invar

####################################

import numpy as np
from dCTES_Strain_functions import *
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


# Define parameters
# Strain of sample (percent)
strain_mech_sample_percent = -0.6

# Calculate sample strain
# Convert strain percentage to strain
strain_mech_sample = strain_mech_sample_percent / 100

# Length of bar (m)
L_sample = 0.1

# Cross-sectional area of bar (m^2)
A_sample = 4.84e-4
A_invar = 2.5e-3
A_aluminum = 1e-2
A_tungsten = 2.5e-3

# Coefficient of thermal expansion (m/m/K)
alpha_sample = 1.2e-5
alpha_invar = 1.85e-6
alpha_aluminum = 2.28e-5
alpha_tungsten = 5.0e-6 #m/m/K https://www.carbideprobes.com/wp-content/uploads/2019/07/TungstenCarbideDataSheet.pdf

# Elastic Modulus (Pa)
E_sample = 130e9
E_invar = 142e9
E_aluminum = 69e9
E_tungsten = 550e9

# Temperature (K)
T_initial = 277.0
T_final = 77.0
dT = T_final - T_initial


#%% Make Graph of sample size vs size of Aluminum and Invar
import matplotlib.pyplot as plt

sample_target_strain = 0.6
deformation_f = -deformation_force_dl(A_sample, L_sample, E_sample, mechanical_deformation_e(L_sample, sample_target_strain))

sample_size = np.arange(0.05, .11,0.01)

all_sizes = np.zeros([len(sample_size),3])
Invar_sizes = all_sizes
for i in range(0,len(all_sizes)):
    all_sizes[i,0] = sample_size[i]
    all_sizes[i,1] = required_length_al(sample_size[i], (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),    
                                        (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
                                        (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)))
    all_sizes[i,2] = all_sizes[i,1]-sample_size[i]
    
fig = plt.figure(1)
ax = fig.gca()
plt.rcParams['figure.dpi'] = 1000

fig.suptitle("Length of Al and In vs Length of Sample for 0.6 % Mechanical Strain",fontsize=10)

#Formatting the axes
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.set_axisbelow(True)

ax.set_xlabel("Sample size [mm]", fontsize=10)
ax.set_ylabel("Aluminum/Invar size [mm]", fontsize=10)

ax.plot(1000*all_sizes[:,0],1000*all_sizes[:,1], label = "Length of Aluminum [mm]")
ax.plot(1000*all_sizes[:,0],1000*all_sizes[:,2], label = "Length of Invar [mm]")

ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='major', color='#CCCCCC', linestyle='--')
ax.grid(which='minor', color='#CCCCCC', linestyle=':')
ax.legend(fontsize = 10)
#%%
L_sample = 0.07
al_trial1 = required_length_al(L_sample, (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),    
                                    (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
                                    (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)))

in_trial1 = al_trial1-L_sample

#%% Strain given a fixed aluminum size

#L_aluminum = al_trial1
L_aluminum = 0.275
L_invar = L_aluminum-L_sample
print(L_invar)

e_th_al = thermal_deformation(alpha_aluminum, L_aluminum, T_initial, T_final)
e_th_s = thermal_deformation(alpha_sample, L_sample, T_initial, T_final)
e_th_i = thermal_deformation(alpha_invar, L_invar, T_initial, T_final)

force_2 = force_from_dL_thermal((e_th_al-e_th_s-e_th_i), L_sample, L_aluminum, L_invar, A_sample , A_aluminum, A_invar, E_sample, E_aluminum, E_invar)
sample_mech_e = mechanical_strain_f(force_2, A_sample, E_sample)

print(deformation_f)
print(force_2)
print("Mechanical strain" + str(round(100*sample_mech_e,4))+"%")


#%%

sample_sizes =np.arange(0.04, 0.11, 0.005)
al_size = 0.1891
all_sizes = np.zeros([len(sample_sizes),4])
#Invar_sizes = all_sizes
fig = plt.figure(1)
#ax1 = fig.gca()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis


for j in range(0,len(all_sizes)):
    all_sizes[j,0] = sample_sizes[j]
    all_sizes[j,1] = al_size
    all_sizes[j,2] = al_size-sample_sizes[j]    
    e_th_al = thermal_deformation(alpha_aluminum, al_size, T_initial, T_final)
    e_th_s = thermal_deformation(alpha_sample, sample_sizes[j], T_initial, T_final)
    e_th_i = thermal_deformation(alpha_invar, al_size-sample_sizes[j], T_initial, T_final)
    force = force_from_dL_thermal((e_th_al-e_th_s-e_th_i), sample_sizes[j], al_size, al_size-sample_sizes[j], A_sample , A_aluminum, A_invar, E_sample, E_aluminum, E_invar)
    sample_mech_e = mechanical_strain_f(force, A_sample, E_sample)
    all_sizes[j,3] = 100*sample_mech_e
    

color = 'tab:red'
ax1.set_xlabel("Sample size [mm]")    
lls = ax1.plot(1000*all_sizes[:,0],all_sizes[:,3], label = str(round(sample_sizes[j],2))+ " mm sample", color = color)#, alpha = 1/(1+j/5))

color = 'tab:blue'
lls2 = ax2.plot(1000*all_sizes[:,0],1000*all_sizes[:,2], label = " Invar size for " + str(round(sample_sizes[j],2)) + " [mm] sample", color = color)#,  alpha = 1/(1+j/5))
    

plt.rcParams['figure.dpi'] = 1000

fig.suptitle("Sample Mechanical Strain vs Sample Size given one Aluminum frame",fontsize=10 )

    #Formatting the axes
fig.tight_layout()
ax1.tick_params(axis='x', labelsize=10)

ax1.tick_params(axis='y', labelsize=10, labelcolor="tab:red")
ax2.tick_params(axis='y', labelsize=10, labelcolor="tab:blue")

ax1.set_ylabel('Sample Mechanical strain %', color="tab:red")
ax2.set_ylabel('Invar size [mm]', color= "tab:blue")  # we already handled the x-label with ax1


ax1.xaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_major_locator(MultipleLocator(10))
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.grid(which='major', color='#CCCCCC', linestyle='--')
ax1.grid(which='minor', color='#CCCCCC', linestyle=':')
#ax2.grid(which='major', color='#CCCCCC', linestyle='solid')
#ax2.grid(which='minor', color='#CCCCCC', linestyle=':')
ax1.legend(fontsize = 8)#, loc = 9)
ax2.legend(fontsize = 8)

#%%
sample_sizes = al_size = np.arange(0.05, 0.11, 0.01)
al_size = np.arange(0.15, .75,0.01)
all_sizes = np.zeros([len(al_size),4])
Invar_sizes = all_sizes
fig = plt.figure(1)
#ax1 = fig.gca()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

for j in range(0,len(sample_sizes)):
    for i in range(0,len(all_sizes)):
        all_sizes[i,0] = sample_sizes[j]
        all_sizes[i,1] = al_size[i]
        all_sizes[i,2] = al_size[i]-sample_sizes[j]
        e_th_al = thermal_deformation(alpha_aluminum, al_size[i], T_initial, T_final)
        e_th_s = thermal_deformation(alpha_sample, sample_sizes[j], T_initial, T_final)
        e_th_i = thermal_deformation(alpha_invar, al_size[i]-sample_sizes[j], T_initial, T_final)
        force = force_from_dL_thermal((e_th_al-e_th_s-e_th_i), sample_sizes[j], al_size[i], al_size[i]-sample_sizes[j], A_sample , A_aluminum, A_invar, E_sample, E_aluminum, E_invar)
        sample_mech_e = mechanical_strain_f(force, A_sample, E_sample)
        all_sizes[i,3] = 100*sample_mech_e
        
    
    color = 'tab:red'
    ax1.set_xlabel("Aluminum size [mm]")    
    lls = ax1.plot(1000*all_sizes[:,1],all_sizes[:,3], label = str(round(sample_sizes[j],2))+ " mm sample", color = color, alpha = 1/(1+j/5))

        
    #ax1 = fig.gca()
    
    color = 'tab:blue'
    lls2 = ax2.plot(1000*all_sizes[:,1],1000*all_sizes[:,2], label = " Invar size for " + str(round(sample_sizes[j],2)) + " [mm] sample", color = color,  alpha = 1/(1+j/5))
        
    
    plt.rcParams['figure.dpi'] = 1000

    fig.suptitle("Sample Mechanical strain vs Aluminum size",fontsize=15)

    #Formatting the axes
fig.tight_layout()
ax1.tick_params(axis='x', labelsize=10)

ax1.tick_params(axis='y', labelsize=10, labelcolor="tab:red")
ax2.tick_params(axis='y', labelsize=10, labelcolor="tab:blue")

ax1.set_ylabel('Sample Mechanical strain %', color="tab:red")
ax2.set_ylabel('Invar size [mm]', color= "tab:blue")  # we already handled the x-label with ax1


ax1.xaxis.set_major_locator(MultipleLocator(50))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_major_locator(MultipleLocator(100))
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
#ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.grid(which='major', color='#CCCCCC', linestyle='--')
ax1.grid(which='minor', color='#CCCCCC', linestyle=':')
#ax2.grid(which='major', color='#CCCCCC', linestyle='solid')
#ax2.grid(which='minor', color='#CCCCCC', linestyle=':')
ax1.legend(fontsize = 5, loc = 9)
ax2.legend(fontsize = 5)


#%%
sample_target_strain = 0.6
deformation_f = -deformation_force_dl(A_sample, L_sample, E_sample, mechanical_deformation_e(L_sample, sample_target_strain))

#sample_size = np.arange(0.05, .26,0.01)

#all_sizes = np.zeros([len(sample_size),3])
#Invar_sizes = all_sizes
#for i in range(0,len(all_sizes)):
#    all_sizes[i,0] = sample_size[i]
sample_size = 0.1 #m
invar_size = 0.170
#aa = required_length_al(sample_size, (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),    
#                                        (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
#                                        (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)))

aa = required_length_al_2(sample_size, invar_size,
                          (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),
                          (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_tungsten, E_tungsten) + thermal_strain(alpha_tungsten, T_initial, T_final)))

print("Al size = " + str(round(aa,3)) + ", In size =" + str(invar_size) + " , Sample size = " + str(sample_size) + ", Tungsten size = " + str(round(aa-(invar_size+sample_size),2)))


#%%
deformation_f = -deformation_force_dl(A_sample, L_sample, E_sample, mechanical_deformation_e(L_sample, sample_target_strain))

#sample_size = np.arange(0.05, .26,0.01)

#all_sizes = np.zeros([len(sample_size),3])
#Invar_sizes = all_sizes
#for i in range(0,len(all_sizes)):
#    all_sizes[i,0] = sample_size[i]
sample_size = 0.1 #m
tungsten_size = 0.02

aa = required_length_al_3(sample_size, tungsten_size,
                          (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),
                          (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_tungsten, E_tungsten) + thermal_strain(alpha_tungsten, T_initial, T_final)))




print("Al size = " + str(round(aa,3)) + ", Tn size =" + str(tungsten_size) + " , Sample size = " + str(sample_size) + ", Invar size = " + str(round(aa-(tungsten_size+sample_size),3)))

#%%

sample_size = np.arange(0.04,0.1,0.01)
tungsten_size = 0.02
results_matirx = np.zeros([len(sample_size),5])

for i in range(0,len(sample_size)):
    results_matirx[i,0] = tungsten_size
    results_matirx[i,1] = sample_size[i]
    deformation_f = -deformation_force_dl(A_sample, sample_size[i], E_sample, mechanical_deformation_e(sample_size[i], sample_target_strain))
    al_size = required_length_al_3(sample_size[i], tungsten_size,
                          (mechanical_strain_f(deformation_f, A_sample, E_sample) + thermal_strain(alpha_sample, T_initial, T_final)),
                          (mechanical_strain_f(-deformation_f, A_aluminum, E_aluminum) + thermal_strain(alpha_aluminum, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_invar, E_invar) + thermal_strain(alpha_invar, T_initial, T_final)), 
                          (mechanical_strain_f(deformation_f, A_tungsten, E_tungsten) + thermal_strain(alpha_tungsten, T_initial, T_final)))
    results_matirx[i,2] = al_size
    results_matirx[i,3] = al_size - tungsten_size - sample_size[i]
    sample_mech_e = mechanical_strain_f(deformation_f, A_sample, E_sample)
    results_matirx[i,4] = sample_mech_e
fig = plt.figure(1)
#ax1 = fig.gca()
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'
ax1.set_xlabel("Sample Size [mm]")    
ax1.plot(1000*results_matirx[:,1],1000*results_matirx[:,2], label = "Aluminum")#, color = color)#, alpha = 1/(1+j/5))
ax1.plot(1000*results_matirx[:,1],1000*results_matirx[:,3], label = "Invar")#, color = color)#, alpha = 1/(1+j/5))
ax1.plot(1000*results_matirx[:,1],1000*results_matirx[:,0], label = "Tungsten Carbide")#, color = color)#, alpha = 1/(1+j/5))

color = 'tab:red'
ax2.plot(1000*results_matirx[:,1],100*results_matirx[:,4], label = "Sample mechanical strain %", color = color)#,  alpha = 1/(1+j/5))
    

plt.rcParams['figure.dpi'] = 1000
fig.suptitle("Required Aluminum and Invar size for different Sample sizes.",fontsize=10 )

    #Formatting the axes
fig.tight_layout()
ax1.tick_params(axis='x', labelsize=10)

ax1.tick_params(axis='y', labelsize=10)#, labelcolor="tab:red")
ax2.tick_params(axis='y', labelsize=10)#, labelcolor="tab:blue")

ax1.set_ylabel('Reqiured sizes [mm]')#, color="tab:red")
ax2.set_ylabel('Sample mechanical strain %')#, color= "tab:blue")  # we already handled the x-label with ax1


ax1.xaxis.set_major_locator(MultipleLocator(10))
ax1.yaxis.set_major_locator(MultipleLocator(50))
ax2.yaxis.set_major_locator(MultipleLocator(0.01))
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
#ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.grid(which='major', color='#CCCCCC', linestyle='--')
ax1.grid(which='minor', color='#CCCCCC', linestyle=':', alpha = 0.5)
#ax2.grid(which='major', color='#CCCCCC', linestyle='solid')
#ax2.grid(which='minor', color='#CCCCCC', linestyle=':')
ax1.legend(fontsize = 8, loc = 1)
ax2.legend(fontsize = 8, loc = 2)

#%%
alum_size = 0.1891
tungsten_size = 0.2
sample_size = 0.07
invar_size = alum_size - tungsten_size - sample_size


e_th_al = thermal_deformation(alpha_aluminum, alum_size, T_initial, T_final)
e_th_s = thermal_deformation(alpha_sample, sample_size, T_initial, T_final)
e_th_i = thermal_deformation(alpha_invar, alum_size-(tungsten_size+sample_size), T_initial, T_final)
e_th_t = thermal_deformation(alpha_tungsten, tungsten_size, T_initial, T_final)
force = force_from_dL_thermal2((e_th_al-e_th_s-e_th_i-e_th_t), sample_size, alum_size, alum_size-(tungsten_size+sample_size), tungsten_size, A_sample, A_aluminum, A_invar, A_tungsten, E_sample, E_aluminum, E_invar, E_tungsten)

sample_mech_e = mechanical_strain_f(force, A_sample, E_sample)
print(sample_mech_e)
#%%

# Calculate the thermal deformation of the sample
dL_thermal_sample = alpha_sample * L_sample * dT

# Calculate the mechanical deformation of the sample
dL_mech_sample = strain_mech_sample * L_sample

# Calculate the total deformation of the sample
dL_total_sample = dL_thermal_sample + dL_mech_sample

# Calculate the force through the assembly using the mechanical deformation of the sample
# Mechanical Deformation = (Force * Length) / (Cross-sectional Area * Elastic Modulus)
# Mechanical Deformation * (Cross-sectional Area * Elastic Modulus) / Length
F = dL_mech_sample * (A_sample * E_sample) / L_sample

# Calculate the strain for the materials
# Thermal Strain = alpha * dT
strain_thermal_sample = alpha_sample * dT
strain_thermal_aluminum = alpha_aluminum * dT
strain_thermal_invar = alpha_invar * dT

# Mechanical Strain = (Force * Length) / (Cross-sectional Area * Elastic Modulus)
# Defined at the start: strain_mech_sample = F / (A_sample * E_sample)
strain_mech_aluminum = -F / (A_aluminum * E_aluminum)
strain_mech_invar = F / (A_invar * E_invar)

# total strain = thermal strain + mechanical strain
strain_total_sample = strain_thermal_sample + strain_mech_sample
strain_total_aluminum = strain_thermal_aluminum + strain_mech_aluminum
strain_total_invar = strain_thermal_invar + strain_mech_invar

# Compatability equations to solve for length of bars
# Deformation of Aluminum = Deformation of Sample + Deformaton of Invar
## 1. L_aluminum * strain_total_aluminum = (L_sample * strain_total_sample) + (L_invar * strain_total_invar)
# Length of Aluminum = Length of Sample + Length of Invar
## 2. L_aluminum = L_sample + L_invar --> L_invar = L_aluminum - L_sample
# Combining Eq. 1 and Eq. 2 to solve the system for L_aluminum:
L_aluminum = (L_sample * (strain_total_sample - strain_total_invar)) / (strain_total_aluminum - strain_total_invar)

# Solve for L_invar by substituting L_aluminum into Eq. 2
L_invar = L_aluminum - L_sample

# Calculate thermal deformations for the other materials
# Thermal deformation = alpha * original length * change in temperature
dL_thermal_aluminum = alpha_aluminum * L_aluminum * dT
dL_thermal_invar = alpha_invar * L_invar * dT

# Calculate mechanical deformations for the other materials
# Mechanical deformation = force * original length / (cross-sectional area * Young's modulus)
# Note: The force on the aluminum is in the opposite direction
dL_mech_aluminum = - F * L_aluminum / (A_aluminum * E_aluminum)
dL_mech_invar = F * L_invar / (A_invar * E_invar)

# Calculate total deformations for the other materials by adding thermal and mechanical deformations
dL_total_aluminum = dL_thermal_aluminum + dL_mech_aluminum
dL_total_invar = dL_thermal_invar + dL_mech_invar

# Print results
print(f"Force: {F:.8f}")

print("\nThermal Strain:")
print(f"Sample: {strain_thermal_sample * 100:.8f}%")
print(f"Aluminum: {strain_thermal_aluminum * 100:.8f}%")
print(f"Invar: {strain_thermal_invar * 100:.8f}%")

print("\nMechanical Strain:")
print(f"Sample: {strain_mech_sample * 100:.8f}%")
print(f"Aluminum: {strain_mech_aluminum * 100:.8f}%")
print(f"Invar: {strain_mech_invar * 100:.8f}%")

print("\nTotal Strain:")
print(f"Sample: {strain_total_sample * 100:.8f}%")
print(f"Aluminum: {strain_total_aluminum * 100:.8f}%")
print(f"Invar: {strain_total_invar * 100:.8f}%")

print("\nTotal Deformation:")
print(f"Sample: {dL_total_sample:.8f}")
print(f"Aluminum: {dL_total_aluminum:.8f}")
print(f"Invar: {dL_total_invar:.8f}")

print("\nLength:")
print(f"Sample: {L_sample:.8f}")
print(f"Aluminum: {L_aluminum:.8f}")
print(f"Invar: {L_invar:.8f}")

print("\nError:")
print(f"dL Error: {(dL_total_aluminum) - (dL_total_invar + dL_total_sample):.8f}")
print(f"L Error: {(L_aluminum) - (L_invar + L_sample):.8f}")

