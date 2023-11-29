import numpy as np

# Define real constants
A_SMALL_NUMBER = 1.0e-10
TWO_PI = 6.283185307
PI = 3.141592654

# Define real arrays for MSU frequencies
MSU_FREQUENCY = np.array([50.299178, 53.740796, 54.960951, 57.949882])  # MSU frequency values
MSU_FREQ = np.array([50.299178, 53.740796, 54.960951, 57.949882])  # Another representation of MSU frequency values
MSU_BANDWIDTH = np.array([0.200, 0.200, 0.200, 0.200])  # MSU bandwidth values
MSU_BANDWIDTHS = np.array([0.200, 0.200, 0.200, 0.200])  # Another representation of MSU bandwidth values
MSU_FREQ_SPLIT_1 = np.array([0.0, 0.0, 0.0, 0.0])  # Frequency split 1 values
MSU_FREQ_SPLIT_2 = np.array([0.0, 0.0, 0.0, 0.0])  # Frequency split 2 values
MSU_EIA = np.array([0.0, 10.71, 21.51, 32.51, 43.91, 56.19])  # MSU EIA values
MSU_NOM_EIAS = np.array([0.0, 10.71, 21.51, 32.51, 43.91, 56.19])  # Nominal EIA values
MSU_VIEW = np.array([0.0, 9.47, 18.94, 28.41, 37.88, 47.35])  # MSU view angles
MSU_VIEW_ANGLES = np.array([0.0, 9.47, 18.94, 28.41, 37.88, 47.35])  # Another representation of MSU view angles
COS2_MSU_VIEW = np.array([1.0, 0.9729, 0.8946, 0.7736, 0.6230, 0.4590])  # COS^2 of MSU view angles

# Define integer array for polarization
MSU_POLARIZATION = np.array([1, 2, 1, 2])  # Polarization values

# Define real constants for scan and orbit time
MSU_SCAN_TIME = 2.962963e-4  # Typical interscan time
MSU_ORBIT_TIME = 0.070799

# Define integer constants
NUM_CHANNELS = 4
NUM_FOVS = 11
N_FOV = 6
NUM_NODES = 2
CENTER_FOV = 6
NUM_VIEW_ANGLES = 5
NUM_SCAN_POSITIONS = 14

# Constants associated with Tb calculations
T_SPACE = 2.73
SPACE = 12
BLACK_BODY = 13

# Constants for various bounds checks
MAX_REASONABLE_TB = np.array([350.0, 350.0, 350.0, 350.0])
MIN_REASONABLE_TB = np.array([150.0, 150.0, 150.0, 150.0])
MAX_REASONABLE_NEAR_NEIGHBOR_DIFF = np.array([15.0, 15.0, 1.5, 1.5])
MAX_HEIGHT = 910.0  # Maximum height of satellite allowed
MIN_HEIGHT = 740.0  # Minimum height of satellite allowed

# Flag values for Ascending/Descending flag
ASCENDING = 1
DESCENDING = 2
UNKNOWN = 0

# Flag Values for Correction Level Flag
RAW_TBS = 0
TBS_CORRECTED_FOR_HEIGHT = 1
TBS_CORRECTED_FOR_DIURNAL_EFFECTS = 2
TBS_CORRECTED_FOR_SAT_TEMP = 3
TBS_MODELLED_NCEP = 4
TBS_CORRECTED_TO_NADIR = 5
TBS_MODELLED_NCEP_DIURNAL = 6
TBS_CORR_FOR_HEIGHT_AND_LAT = 7
TBS_CORR_TO_NOM_EIA = 8

# Pseudo-enumerated types
SURF_TYPE_LAND = 1
SURF_TYPE_OCEAN = 0
