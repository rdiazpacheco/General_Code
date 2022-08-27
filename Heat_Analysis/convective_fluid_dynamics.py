"""
Supercurrent wand Helium convection cooling calculations
last modified 1/11/21
author: JLC 
"""

from math import pi
import Constants
import materials_properties

def m_dot():
	""" quick computation for mass flow rate given constants
	Returns:
		m_dot: [float] mass flow rate [g/s]
	"""
	_,density_RT,_,_=materials_properties.properties_He(Constants.ROOM_TEMP_K,Constants.SET_PRESSURE_MOCHI) # density [g/l]
	# print("density at circulation pump = "+str(density_RT) + " g/l")

	m_dot=round(Constants.FLOW_RATE_MOCHI/Constants.MIN_TO_SEC*density_RT,Constants.ROUND_PREC) # [g/s]
	# print("m_dot = "+str(m_dot) + " g/s")

	return m_dot


def hydraulic_diameter(radius):
	"""
	Args:
		radius: [float] length [m]
	Returns:
		Ac: [float] cross sectional area [m^2]
		Dh: [float] hydraulic diameter [m]	
	"""
	Ac= pi*(radius)**2 # cross_sectional_area [m^2]
	per = 2*pi*radius # wetted perimeter for round shape [m]
	Dh=4*Ac /per # hydraulic diameter [m]

	return Ac, Dh


def reynolds_number(t,radius_inner):
	""" Returns reynolds number for probe in insert, approximated as annulus (concentric pipes)
	https://www.sciencedirect.com/topics/earth-and-planetary-sciences/annular-flow#:~:text=According%20to%20a%20publication%20of,pipe%20in%20the%20classic%20calculation.
	# reference equations 5-3, 5-7, 5-8 in Heat Transfer by Gregory Nellis and Sanford Klein
	we can assume that flow is incompressible if the divergence =0, also generally true at slow speeds
	Args:
		t: [float] helium temperature [K]
		radius_probe: [float] length [m]
	Returns:
		[float] unitless. Laminar flow if > Re_critical
	"""
	Ac_tube,Dh_tube= hydraulic_diameter(Constants.INSERT_RADIUS)
	Ac_probe,Dh_probe=hydraulic_diameter(radius_inner)

	cp,rho,mu,_ =materials_properties.properties_He(t,Constants.SET_PRESSURE_MOCHI) # Cp [J/g-K], density [g/l], viscosity [g/(m-s)]
	# print("He Specific heat at point of interest = "+str(cp)+" J/g-K")
	print("He density at point of interest = "+str(rho)+" g/l")
	print("He viscosity = "+str(mu)+" g/(m*s)")
	area=Ac_tube-Ac_probe
	print("Area = "+str(area)+" m^2")
	u_m= m_dot()/(rho*(area)) # bulk velocity [l/(s-m^2)]
	print("bulk velocity ="+str(u_m)+" l/(s*m^2)") 
	u=u_m*0.001
	print("bulk velocity ="+str(u)+" m/(s)")
	dh=Dh_tube-Dh_probe # m
	print("Dh = "+str(dh)+" m")
	reynolds=(rho*dh*u_m)/mu 
	print("reynolds number = "+str(reynolds))

	return reynolds


def prandtl_number(t):
	""" returns prandtl number, which is the kinematic viscosity / thermal diffusivity 
	Pr=(viscosity/density) / (thermal conductivity/(density*specific heat capacity))= mu * Cp / K 
	Args:
		t: [float] helium temperature [K]
	Returns:
		[float] unitless
	"""
	cp,rho,mu,k =materials_properties.properties_He(t,Constants.SET_PRESSURE_MOCHI) # Cp [J/g-K], density [g/l], viscosity [g/(m-s)], k [W/(m-K)]

	return mu*cp/ k


def local_heat_transfer_coefficient(t,radius):
	""" returns local thermal transfer coefficient, which is = Nusselt number * thermal conductivity /hydraulic diameter
	reference equations 5-39, in Heat Transfer by Gregory Nellis and Sanford Klein
	Args:
		t: [float] helium temperature [K]
		radius: [float] radius of probe [m]
	Returns:
		h: [float] local thermal transfer coefficien [W/(m^2-K)]
	"""
	_,_,_,k =materials_properties.properties_He(t,Constants.SET_PRESSURE_MOCHI) # k [W/(m-K)]

	Dh=hydraulic_diameter(Constants.INSERT_RADIUS)[1]-hydraulic_diameter(radius)[1]
	h=Constants.NUSSELT_NUMBER_ANNULUS_LAMINAR*k/Dh
	print("local thermal transfer coefficient"+h)
	return h