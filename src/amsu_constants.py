
import numpy as np

# Define real constants
A_SMALL_NUMBER = 1.0e-10
TWO_PI = 6.283185307
PI = 3.141592654

# Define real arrays for AMSU frequencies
AMSU_A_FREQ = np.array([23.800, 31.400, 50.300, 52.800, 53.596, 54.400, 54.940, 55.500,
                        57.290334, 57.290334, 57.290334, 57.290334, 57.290334, 57.290334, 89.000])  # Channel frequencies

AMSU_A_STOPBAND = np.array([0.018, 0.018, 0.018, 0.018, 0.0, 0.018, 0.018, 0.018, 0.018,
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # Stopband values

AMSU_A_FREQ_SPLIT_1 = np.array([0.0, 0.0, 0.0, 0.0, 0.115, 0.0, 0.0, 0.0, 0.0, 0.217,
                                0.3222, 0.3222, 0.3222, 0.3222, 0.0])  # Frequency split 1 values

AMSU_A_FREQ_SPLIT_2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.048,
                                0.022, 0.010, 0.0045, 0.0])  # Frequency split 2 values

AMSU_A_BANDWIDTH = np.array([0.251, 0.161, 0.161, 0.3805, 0.170, 0.3805, 0.3805, 0.3103,
                             0.3300, 0.07658, 0.03511, 0.01529, 0.00793, 0.00294, 1.9989])  # Bandwidth values

# Define integer constant
NUM_AMSU_FOVS = 30  # Number of AMSU FOVs

# Define real arrays for view angles and nominal EIAS
AMSU_VIEW_ANGLES = np.array([1.6666666, 5.0000000, 8.3333333, 11.666667, 15.000000,
                             18.333333, 21.666667, 25.000000, 28.333333, 31.666667,
                             35.000000, 38.333333, 41.666667, 45.000000, 48.333333])  # View angles

AMSU_NOM_EIAS = np.array([1.875947, 5.629541, 9.388301, 13.155880, 16.936250,
                          20.733890, 24.554020, 28.402830, 32.287970, 36.219100,
                          40.208800, 44.274040, 48.438570, 52.737200, 57.224260])  # Nominal EIAS values

# Define integer array for polarization
AMSU_A_POLARIZATION = np.array([1, 1, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1])  # Polarization values

AMSU_A_TLT_WTS = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.25,0.40,1.17,1.61,1.41,0.44,-1.14,-2.64])
AMSU_A_TMT_WTS = np.array([1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
AMSU_A_TMT_WTS = AMSU_A_TMT_WTS/np.sum(AMSU_A_TMT_WTS)


x = 1