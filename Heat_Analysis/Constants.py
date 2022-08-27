"""
Constants and Enums
author: JLC 
"""

#%%
import enum 

# physical constants

ROOM_TEMP_K = 298 	# K
TEMP_S1_K = 40 		# K
TEMP_S2_K = 6.8 	# K
TEMP_SAMPLE_K = 20 	# K
TEMP_LOWER_K = 60	# K, range 50-75 K
TEMP_RECUP_K = 80	# K, guess

MAX_Q_S1_W = 7.5 	# W
MAX_Q_S2_W = 16 	# W

FLOW_RATE_MOCHI=4.5 # l/min
SET_PRESSURE_MOCHI= 10 # psi
INSERT_RADIUS = 0.0106 # m
PROBE_RADIUS= 0.01896/2 # m, note this is an approximation
REYNOLDS_CRITICAL=300000

"""
 assumptions for NUSSELT Number Annular duct, p. 667 Heat Transfer by Gregory Nellis and Sanford Klein
adiabatic external surface because insert is insulated, 
internal surface constant heat flux, (steady state?)
Radius ratio = Probe radius/ insert radius = 0.894
"""
NUSSELT_NUMBER_ANNULUS_LAMINAR= 4.8

# Conversion factors

kW_TO_W=1000 # convert from kW to W
kG_TO_G=1000 # convert from kg to g
HE_MOL_TO_G=4.0 # convert from mols to grams
MIN_TO_SEC=60 # convert from minutes to seconds
uPA_TO_G_PER_M_S2= 0.001 # convert from uPa to g/(m*s^2), note 1 Pa = 1 Kg / (m*s^2)
uOHM_CM_TO_OHM_M=1e-8 

# paths
HE_PROPERTIES_PATH_BASE="HE-isobaric-properties-overpressure-psi-"


# other
ROUND_PREC=4

# enums
class MATERIAL(enum.Enum):
    CU_R3_50="OFHC Copper RRR=50"
    CU_R3_100="OFHC Copper RRR=100"
    CU_R3_150="OFHC Copper RRR=150"
    SS_316="316 Stainless"
    G_10="G10"
    HTS="HTS"

class GEOMETRY(enum.Enum):
	ANNULAR="annular"
