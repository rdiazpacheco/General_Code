"""
title
author: JLC 
"""
import heat_analysis
import convective_fluid_dynamics
import materials_properties
import pytest
import Constants

# ########################################Test Materials Properties #################################

def test_thermal_conductivity_of():
	# test values compared to Ekin A3.1 p.517 thermal conductivity 
	assert round(materials_properties.thermal_conductivity_of(Constants.MATERIAL.CU_R3_50.value,293)/1000,2)==round(0.394,2) # kW/(m-K)
	# assert round(materials_properties.thermal_conductivity_of(Constants.MATERIAL.CU_R3_150.value,4.2)/1000,2)==round(0.850,2) # kW/(m-K)

def test_properties_He():
	# test values compared to exported to NIST data folder from :
	# https://webbook.nist.gov/cgi/fluid.cgi?ID=C7440597&TUnit=K&PUnit=atm&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&Type=IsoBar&RefState=DEF&Action=Page
	with pytest.raises(Exception):
		materials_properties.properties_He(300,4.9)
	with pytest.raises(Exception):
		materials_properties.properties_He(300,10.1)
	with pytest.raises(Exception):
		materials_properties.properties_He(300.1,5)
	with pytest.raises(Exception):
		materials_properties.properties_He(3.8,10.5)
	assert materials_properties.properties_He(300,5)[0]==5.19625
	assert materials_properties.properties_He(300,5)[1]==0.21764

	assert materials_properties.properties_He(20.2,10)[0]==5.29325
	assert materials_properties.properties_He(20.2,10)[1]==4.0988

	assert round(materials_properties.properties_He(20,10)[0],2)==5.29
	assert round(materials_properties.properties_He(20,10)[1],2)==4.10
	assert round(materials_properties.properties_He(20,10)[2],4)==0.0036
	assert round(materials_properties.properties_He(20,10)[3],4)==0.0263

def test_specific_heat_of(): # J/g-K
	# test values compared against Ekin A6.2 Specific heat vs. Temperature on page 569
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_50.value,300),2)==round(0.386,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_100.value,300),2)==round(0.386,2) # J/g-K , use same curve
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_150.value,300),2)==round(0.386,2) # J/g-K , use same curve
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_50.value,4),2)==round(0.00009,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_50.value,20),2)==round(0.0070,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.CU_R3_50.value,77),1)==round(0.192,1) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.SS_316.value,300),1)==round(0.48,1) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.SS_316.value,4),2)==round(0.0020,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.SS_316.value,20),1)==round(0.017,1) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.SS_316.value,77),2)==round(0.20,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.G_10.value,300),2)==round(0.999,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.G_10.value,4),2)==round(0.0020,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.G_10.value,20),2)==round(0.047,2) 
	assert round(materials_properties.specific_heat_of(Constants.MATERIAL.G_10.value,77),2)==round(0.239,2) 

def test_resistivity_of(): # Ohm-m
	# Data from Appendix 6.5a and b in Ekin p.575-576, in Ohm*m*10^-8   
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,10),8)==round(53.9*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,20),8)==round(53.9*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,50),8)==round(54.9*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,77),8)==round(56.8*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,100),8)==round(58.8*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,200),8)==round(68.9*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.SS_316.value,295),8)==round(77.1*1e-8,8) 

	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,10),8)==round(0.015*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,20),8)==round(0.017*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,50),8)==round(0.084*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,77),8)==round(0.21*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,100),8)==round(0.34*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,200),8)==round(1.07*1e-8,8) 
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.CU_R3_100.value,295),8)==round(1.70*1e-8,8) 		

	assert round(materials_properties.resistivity_of(Constants.MATERIAL.HTS.value,20),2)==round(0,2)
	assert round(materials_properties.resistivity_of(Constants.MATERIAL.HTS.value,100),2)==round(0.3957*1e-8,2)


# ########################################### Test Convective Fluid Analysis ################################################

def test_m_dot():
	# cross referenced against "A 1 kA-class cryogen-free critical current characterization system for superconducting coated conductors" Strickland
	# under section IV model of gas flow cooling = 0.03 g/s at 7.8 L/min
	assert convective_fluid_dynamics.m_dot()==0.0206

def test_reynolds_number():
	# check that reynolds number is calculated correctly throughout temperature range
	assert round(convective_fluid_dynamics.reynolds_number(12,0.001),2)==434.39
	assert round(convective_fluid_dynamics.reynolds_number(20,Constants.PROBE_RADIUS),2)==180.94
	assert round(convective_fluid_dynamics.reynolds_number(88,Constants.PROBE_RADIUS),2)==72.30
	assert round(convective_fluid_dynamics.reynolds_number(275,Constants.PROBE_RADIUS),2)==34.77


# ########################################### Test Heat Analysis ###########################################################

def test_thermal_conductivity_integral():
	# Data from Appendix A2.1 in Ekin p. 514-515 in kW/m

	# approximate ETP Cu as RRR=50, divide numbers to be in W/m with multiplied factors for roundoff to 2 sig figs 
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_50.value,4,20)/1e5,2)==round(14*1e-2,2) 
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_50.value,4,70)/1e5,1)==round(65.1*1e-2,1)
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.CU_R3_50.value,4,300)/1e5,1)==round(162*1e-2,1) 
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.SS_316.value,4,20)/1e3,1)==round(0.0163,1) 
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.SS_316.value,4,70)/1e3,1)==round(0.270,1)
	assert round(heat_analysis.thermal_conductivity_integral(Constants.MATERIAL.SS_316.value,4,300)/1e4,1)==round(3.06*1e-1,1) 


def test_heat_transfer_rate():
	assert heat_analysis.heat_transfer_rate(Constants.TEMP_S2_K,Constants.TEMP_S1_K,Constants.FLOW_RATE_MOCHI,Constants.SET_PRESSURE_MOCHI)>=-16


