####################################
# Thermal Stress Test V2
# Solve for L_aluminum and L_invar
####################################
import numpy as np
import matplotlib.pyplot as plt
# Define parameters
# Strain of sample (percent)
strain_sample_percent = -0.6

# Calculate sample strain
# Convert strain percentage to strain
strain_sample = strain_sample_percent / 100

# Length of bar (m)
L_sample = 0.05
L_sample_array = np.linspace(0.015, 0.05,100)

# Cross-sectional area of bar (m^2)
A_sample = 4.84e-4
A_invar = 2.5e-3
A_aluminum = 1e-2

# Coefficient of thermal expansion (1/K)
alpha_sample = 1.2e-5
alpha_invar = 1.85e-6
alpha_aluminum = 2.28e-5

# Elastic Modulus (Pa)
E_sample = 130e9
E_invar = 142e9
E_aluminum = 69e9

# Temperature (K)
T_initial = 277.0
T_final = 77.0
dT = abs(T_final - T_initial)

# Calculate the total deformation of the sample using strain
dL_total_sample = strain_sample * L_sample
    
# Calculate the thermal deformation of the sample
dL_thermal_sample = - alpha_sample * L_sample * dT

# Calculate the mechanical deformation of the sample
dL_mech_sample = dL_total_sample - dL_thermal_sample
    
# Calculate the force using the mechanical deformation of the sample
# Mechanical deformation = FL/AE
F = dL_mech_sample * (A_sample * E_sample) / L_sample

# Compatability equations to solve for length of bars
# Deformation of Aluminum = Deformation of Sample + Deformaton of Invar
## 1. L_aluminum*(-F/(A_aluminum*E_aluminum)+(alpha_aluminum*dT)) = L_sample*(F/(A_sample*E_sample)+(alpha_sample*dT)) + L_invar*(F/(A_invar*E_invar)+(alpha_invar*dT))
# Length of Aluminum = Length of Sample + Length of Invar
## 2. L_aluminum = L_sample + L_invar --> L_invar = L_aluminum - L_sample
# Combining Eq. 1 and Eq. 2 to solve the system for L_aluminum:
L_aluminum = (A_aluminum * E_aluminum * (A_invar * alpha_invar * A_sample * dT * E_invar * E_sample - A_invar * alpha_sample * A_sample * dT * E_invar * E_sample + A_invar * E_invar * F - A_sample * E_sample * F) * L_sample) / (A_sample * E_sample * (-A_aluminum * A_invar * alpha_aluminum * dT * E_aluminum * E_invar + A_aluminum * A_invar * alpha_invar * dT * E_aluminum * E_invar - A_aluminum * E_aluminum * F - A_invar * E_invar * F))


#Behavior - 
L_aluminum_array = (A_aluminum * E_aluminum * (A_invar * alpha_invar * A_sample * dT * E_invar * E_sample - A_invar * alpha_sample * A_sample * dT * E_invar * E_sample + A_invar * E_invar * F - A_sample * E_sample * F) * L_sample_array) / (A_sample * E_sample * (-A_aluminum * A_invar * alpha_aluminum * dT * E_aluminum * E_invar + A_aluminum * A_invar * alpha_invar * dT * E_aluminum * E_invar - A_aluminum * E_aluminum * F - A_invar * E_invar * F))

L_aluminum2 = (L_sample*((dT*(alpha_sample-alpha_invar))+(F*((1/(A_sample*E_sample))-(1/(A_invar*E_invar))))))/((dT*(alpha_aluminum-alpha_invar))-(F*((1/(A_aluminum*E_aluminum))+(1/(A_invar*E_invar)))))

#L_aluminum2 = (L_sample*((dT*(-alpha_sample+alpha_invar))+(F*((1/(A_sample*E_sample))-(1/(A_invar*E_invar))))))/((dT*(-alpha_aluminum+alpha_invar))-(F*((1/(A_aluminum*E_aluminum))+(1/(A_invar*E_invar))))) #this one matches Steven's

M = A_aluminum*E_aluminum
N = A_invar*E_invar
O = A_sample*E_sample

L_aluminum3 = (L_sample*M*((F*(O-N))+(N*O*dT*(alpha_invar-alpha_sample))))/(O*((F*(M+N))+(M*N*dT*(alpha_invar-alpha_aluminum))))


plt.plot(L_sample_array,L_aluminum_array)


# Solve for L_invar by substituting L_aluminum into Eq. 2
L_invar = L_aluminum - L_sample

# Calculate thermal deformations for the other materials
# Thermal deformation = alpha * original length * change in temperature
dL_thermal_invar = - alpha_invar * L_invar * dT
dL_thermal_aluminum = - alpha_aluminum * L_aluminum * dT

# Calculate mechanical deformations for the other materials
# Mechanical deformation = force * original length / (cross-sectional area * Young's modulus)
# Note: The force on the aluminum is in the opposite direction
dL_mech_invar = F * L_invar / (A_invar * E_invar)
dL_mech_aluminum = - F * L_aluminum / (A_aluminum * E_aluminum)

# Calculate total deformations for the other materials by adding thermal and mechanical deformations
dL_total_invar = dL_thermal_invar + dL_mech_invar
dL_total_aluminum = dL_thermal_aluminum + dL_mech_aluminum

# Calculate strain for the other materials
# Strain = total deformation / original length * 100%
strain_invar_percent = (dL_total_invar / L_invar) * 100
strain_aluminum_percent = (dL_total_aluminum / L_aluminum) * 100

# Print results
print(f"Force: {F:.8f}")

print("\nThermal Deformation:")
print(f"Sample: {dL_thermal_sample:.8f}")
print(f"Aluminum: {dL_thermal_aluminum:.8f}")
print(f"Invar: {dL_thermal_invar:.8f}")

print("\nMechanical Deformation:")
print(f"Sample: {dL_mech_sample:.8f}")
print(f"Aluminum: {dL_mech_aluminum:.8f}")
print(f"Invar: {dL_mech_invar:.8f}")

print("\nTotal Deformation:")
print(f"Sample: {dL_total_sample:.8f}")
print(f"Aluminum: {dL_total_aluminum:.8f}")
print(f"Invar: {dL_total_invar:.8f}")

print("\nStrain:")
print(f"Sample: {strain_sample_percent:.8f}%")
print(f"Aluminum: {strain_aluminum_percent:.8f}%")
print(f"Invar: {strain_invar_percent:.8f}%")

print("\nError:")
print(f"L Error: {(L_aluminum) - (L_invar + L_sample):.8f}")
print(f"dL Error: {(dL_total_aluminum) - (dL_total_invar + dL_total_sample):.8f}")
