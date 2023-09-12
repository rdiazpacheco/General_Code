####################################
# Thermal Stress Test V2
# Solve for L_aluminum and L_invar
####################################

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
dT = T_final - T_initial
    
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