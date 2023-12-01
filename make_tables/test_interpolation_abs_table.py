#program test_interpolation_abs_table
import datetime
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import xarray as xr
import time

import atm_rtm
import test.msu_constants as msu_constants
import test.amsu_constants as amsu_constants
import matplotlib.pyplot as plt

from amsu_rtm import AMSU_O2_Q_INTERPOLATOR, AMSU_CLD_INTERPOLATOR

channel = 5
AMSU_num_freq = 7
O2_model = 'RSS_2022'
o2_q_interp = AMSU_O2_Q_INTERPOLATOR(channel, O2_model)
cld_interp = AMSU_CLD_INTERPOLATOR(channel)

p_list  = np.arange(10.0,1010.0,0.1,dtype=np.float32)
t_list = np.arange(200.0,320.0,0.1,dtype=np.float32)
abs1_array=np.zeros_like(p_list)
abs2_array=np.zeros_like(p_list)

interp_points = np.zeros((3,len(p_list)),dtype=np.float32,order='F')
abs_at_interp_points = np.zeros((len(p_list)),dtype=np.float32,order='F')
for ip,p_value in enumerate(p_list):
    interp_points[:,ip] = np.array([0.010001, p_value, 250.0001])


# calculate the absorption coefficient for interp_points using the interpolator
st_cpu = time.process_time()
st = time.time()
for i in range(10):
    print(i)
    abs1_array = o2_q_interp(interp_points[0,:],interp_points[1,:],interp_points[2,:])

et_cpu = time.process_time()
et = time.time()
res = et-st
print(f'time:RegularGridInterpolator:amsu_rtm:AMSU_O2_Q_INTERPOLATOR {res}')
res = et_cpu-st_cpu
print(f'cpu_time:RegularGridInterpolator:amsu_rtm:AMSU_O2_Q_INTERPOLATOR {res}')


    
# calculate the absorption coefficient for a single point directly
AMSU_num_freq = 7
channel_index = channel - 1
if amsu_constants.AMSU_A_FREQ_SPLIT_2[channel_index] > 0.01:
    raise ValueError('This code does not support AMSU-A channels with double-split frequency bands')

if amsu_constants.AMSU_A_FREQ_SPLIT_1[channel_index] < 0.01:
    AMSU_num_band = 1
else:
    AMSU_num_band = 2

# construct an array of frequencies where the calc is performed
amsu_freq_arr = np.full((AMSU_num_freq,AMSU_num_band),np.nan,dtype=np.float32,order='F')
for band_index in range(0,AMSU_num_band):
    for freq_index in range(0,AMSU_num_freq):
        amsu_freq_arr[freq_index,band_index] = (                        
                amsu_constants.AMSU_A_FREQ[channel_index] +                                    
                2.0*(band_index-0.5)*amsu_constants.AMSU_A_FREQ_SPLIT_1[channel_index] -                 
                amsu_constants.AMSU_A_BANDWIDTH[channel_index]/2.0 +                            
                (1+2*freq_index)*amsu_constants.AMSU_A_BANDWIDTH[channel_index]/(2.0*AMSU_num_freq))
amsu_freq_arr = amsu_freq_arr.flatten(order='F')
for ip,p_value in enumerate(p_list):
    #print(amsu_freq_arr)
    
    o2abs = np.zeros((1),dtype=np.float32)
    h2oabs = np.zeros((1),dtype=np.float32)

    abs_tot = 0.0
    t_test = interp_points[2,ip]
    p_test = interp_points[1,ip]
    q_test = interp_points[0,ip]
    for freq in amsu_freq_arr:
        pv = p_test*q_test/(0.622+0.378*q_test)
        atm_rtm.atm_rtm.abs_o2_rss_2022(t_test,p_test,pv,freq,o2abs)
        abs_tot = abs_tot + o2abs
        atm_rtm.atm_rtm.abs_h2o_rss_2022(t_test,p_test,pv,freq,h2oabs)
        abs_tot = abs_tot + h2oabs
        #print(f'freq: {freq} o2abs: {o2abs} h2oabs: {h2oabs} abs_tot: {o2abs+h2oabs}')

    abs = abs_tot/len(amsu_freq_arr)
    abs2_array[ip] = abs[0]

# compare the two results

print(abs1_array)
print(abs2_array)
print((abs1_array-abs2_array))

fig1,ax1 = plt.subplots(1,1)
ax1.plot(p_list,abs1_array)
ax1.plot(p_list,abs2_array)
ax1.set_xlabel('Pressure (hPa)')
ax1.set_ylabel('Absorption Coefficient (nepers/km$)')

fig2,ax2 = plt.subplots(1,1)
ax2.patch.set_alpha(1.0)
ax2.patch.set_facecolor('white')
ax2.plot(p_list,(abs1_array-abs2_array)/abs2_array)
ax2.set_ylim(0.0,0.001)
ax2.set_xlabel('Pressure (hPa)')
ax2.set_ylabel('Fractional Difference')


# calculate the cloud absorption coefficient for interp_points using the interpolator

cld_list = np.arange(0.0,0.3,0.001,dtype=np.float32)

abs1_array=np.zeros_like(cld_list)
abs2_array=np.zeros_like(cld_list)

interp_points = np.zeros((2,len(cld_list)),dtype=np.float32,order='F')

for ip,cld_value in enumerate(cld_list):
    interp_points[:,ip] = np.array([cld_value,250.0001])

abs_at_interp_points = cld_interp(interp_points[0,:],interp_points[1,:])
abs1_array = abs_at_interp_points

# calculate the cloud absorption coefficient for a single point directly


for ip,cld_value in enumerate(cld_list):
    for freq in amsu_freq_arr:
        abs_cld = np.zeros((1),dtype=np.float32)
        abs_cld = atm_rtm.atm_rtm.fdcldabs(freq,interp_points[1,ip],interp_points[0,ip])
        abs2_array[ip] += abs_cld
abs2_array = abs2_array/(AMSU_num_freq*AMSU_num_band)

# compare the two results
print(abs1_array)
print(abs2_array)
print((abs1_array-abs2_array))

fig1,ax1 = plt.subplots(1,1)
ax1.plot(cld_list,abs1_array)
ax1.plot(cld_list,abs2_array)
ax1.set_xlabel('Cloud Density (kg/m$^3$)')
ax1.set_ylabel('Absorption Coefficient (nepers/km$)')
ax1.set_title(f'AMSU Channel {channel}')

fig2,ax2 = plt.subplots(1,1)
ax2.plot(cld_list,(abs1_array-abs2_array)/abs2_array)
ax2.set_ylim(0.0,0.001)
ax2.set_xlabel('Cloud Density (kg/m$^3$)')
ax2.set_ylabel('Fractional Difference')
ax2.set_title(f'AMSU Channel {channel}')

print()
