
import numpy as np
import os
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import xarray as xr

#from atm_abs import o2abs_19_sub, absoxy_rss_2022_sub
import atm_rtm

import sys
sys.path.append('/mnt/m/job_access/atm_abs_python_fortran/make_tables')
import test.msu_constants as msu_constants
import test.amsu_constants as amsu_constants

import matplotlib.pyplot as plt
from rss_plotting.global_map import plot_global_map

rd=287.05
epsilon=0.622
rv=rd/epsilon
one_minus_epsilon=1.0-epsilon
epsilon_ratio=one_minus_epsilon/epsilon

def convert_ql_to_rhol(ql,p,q,t):
    # ql: cloud liquid water (kg/kg)
    # p: pressure (hPa)
    # q: specific humidity (kg/kg)
    # t: temperature (K)
    # returns: cloud liquid water density (kg/m^3)

    # based on the treatment in 
    # https://www.nwpsaf.eu/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
    # section 4. Retrieved 11/23/2023


    rm = rd*(1+epsilon_ratio*q)
    # calculate rho cloud liquid water
    rhol = ql*100.0*p/(rm*t)
    return rhol

def calc_tb_toa(tran,tbup,tbdw,t_surf,emiss,space_tb):
    return tbup + tran*(t_surf*emiss + (1.0-emiss)*(tbdw+(tran*space_tb)))

def find_frequencies(channel):
    channel_index = channel - 1
    AMSU_num_freq = 7
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


    return amsu_freq_arr.flatten(order='F')

re = 6371.0e3
amsu_channel = 5
O2_model = 'RSS_2022'
nc_file = Path(f'/mnt/m/job_access/atm_abs_python_fortran/make_tables/tables/amsu_{amsu_channel:02d}_abs_table_q_per_km_{O2_model}.v3.nc')
print(nc_file)
ds = xr.open_dataset(nc_file)
ds.keys()
abs_table = ds['absorptivity'].values
q_table_axis = ds['specific_humidity'].values
p_table_axis = ds['pressure'].values
t_table_axis = ds['temperature'].values

era5_file = Path('/mnt/a/data/_access_temp/rtm/era5_levels_2022-04-16.nc')
era5_surf_file = Path('/mnt/a/data/_access_temp/rtm/era5_surface_2022-04-16.nc')

ds = xr.open_dataset(era5_file)
ds_surf = xr.open_dataset(era5_surf_file)

ilat = 360
ilon = 0
time = 0
p = ds['level'].values.astype(np.float32)
p = np.asfortranarray(np.tile(p,(1440,721,1)).T)
t = np.asfortranarray(ds['t'][time,:,:,:].values)
ql = np.asfortranarray(ds['clwc'][time,:,:,:].values)
q = np.asfortranarray(ds['q'][time,:,:,:].values)
z = np.asfortranarray(ds['z'][time,:,:,:].values)
p_surf = np.asfortranarray(ds_surf['sp'][time,:,:].values/100.0)
t_surf = np.asfortranarray(ds_surf['t2m'][time,:,:].values)
dp_surf = np.asfortranarray(ds_surf['d2m'][time,:,:].values)
z_surf = np.asfortranarray(ds_surf['z'][time,:,:].values)
h_surf =z_surf*re/(9.80665*re - z_surf)

pv_surf = np.asfortranarray(np.zeros((721*1440),dtype=np.float32))

atm_rtm.atm_rtm.goff_gratch_rss_2022(np.reshape(dp_surf,-1),
                                    np.ones((721*1440),dtype=np.float32)*100.0,
                                    np.reshape(p_surf,-1),
                                    pv_surf)
# flatten the lat/lon dimensions
p_surf = np.reshape(p_surf,-1)
t_surf = np.reshape(t_surf,-1)
dp_surf = np.reshape(dp_surf,-1)
h_surf = np.reshape(h_surf,-1)
# pv_surf is already flattened

t = np.reshape(t,(t.shape[0],-1))
p = np.reshape(p,(p.shape[0],-1))
ql = np.reshape(ql,(ql.shape[0],-1))
q = np.reshape(q,(q.shape[0],-1))
z = np.reshape(z,(z.shape[0],-1))

#convert from geopontential to height
h =z*re/(9.80665*re - z)

#convert from q to pv
pv = p*q/(0.622+0.378*q)

# add the surface values to the bottom of the profiles
h = np.insert(h,0,h_surf,axis=0)
p = np.insert(p,0,p_surf,axis=0)
t = np.insert(t,0,t_surf,axis=0)
q = np.insert(q,0,0.0,axis=0)
pv = np.insert(pv,0,pv_surf,axis=0)
ql = np.insert(ql,0,0.0,axis=0)

rhol = convert_ql_to_rhol(ql,p,q,t)

shp = t.shape
nlev = shp[0]
nprofile = shp[1]
# amsu_channel = 5
# center_freq_list = [53.48,53.71]
# bandwidth_list = [0.17,0.17]
# n_sub_band = 5

# find the cloud absorption 

amsu_frequencies = find_frequencies(amsu_channel)

al = np.zeros_like(p,dtype=np.float32,order='F')
al_tot = np.zeros_like(p,dtype=np.float32,order='F')
for freq in amsu_frequencies:
    print(f'Processing freq {freq:.2f} GHz')
    atm_rtm.atm_rtm.fdcldabs_2d_ql(freq,p,t,q,ql,al)
    al_tot = al_tot + al

cld_abs = al_tot/len(amsu_frequencies)

# find the clear sky absorption

interp_method = RegularGridInterpolator((q_table_axis,p_table_axis,t_table_axis),abs_table,
                                        method='linear',bounds_error=False,fill_value=None)
interp_points = np.array([q,p,t]).T
print('done with init for RegularGridInterpolator')
st_cpu = time.process_time()
st = time.time()
abs1_array = interp_method(interp_points)
et_cpu = time.process_time()
et = time.time()
res = et-st
print(f'time:RegularGridInterpolator: {res}')
res = et_cpu-st_cpu
print(f'cpu_time:RegularGridInterpolator: {res}')

'''THE ABOVE IS NOT TESTED YET  
ALSO NEED TO ADD THE CLOUD ABSORPTION TO THE CLEAR SKY ABSORPTION
AND CALL atm_tran_2d'''
    


tb_maps_rss = np.zeros((2,5,721,1440),dtype=np.float32,order='F')
tb_maps_rosen = np.zeros((2,5,721,1440),dtype=np.float32,order='F')

for side_band,center_freq in enumerate(center_freq_list):
    delta_freq = bandwidth_list[side_band]/n_sub_band
    k_max = np.floor(n_sub_band/2)
    freq_list = np.arange(-k_max,k_max+1)*delta_freq + center_freq
    for ifreq,freq in enumerate(freq_list):
        print(f'Processing freq {freq:.2f} GHz')
        theta = 0.0
        emiss = 0.9
        space_tb = 2.73

        num_profiles = h.shape[1]

        tran_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
        tbup_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
        tbdw_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')

        tran_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
        tbup_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
        tbdw_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')

        atm_rtm.atm_rtm.get_atm_components_rss_2022_2d(
                    freq, 
                    theta, 
                    h, 
                    p, 
                    t, 
                    pv, 
                    clwc,
                    tran_rss,
                    tbup_rss,
                    tbdw_rss)

        atm_rtm.atm_rtm.get_atm_components_rosen_2017_2d(
                    freq, 
                    theta, 
                    h, 
                    p, 
                    t, 
                    pv, 
                    clwc,
                    tran_rosen,
                    tbup_rosen,
                    tbdw_rosen)



        tb_toa_rss = calc_tb_toa(tran_rss,tbup_rss,tbdw_rss,t_surf,emiss,space_tb)
        tb_toa_rss = np.reshape(tb_toa_rss,(721,1440))
        tb_maps_rss[side_band,ifreq,:,:] = tb_toa_rss

        tb_toa_rosen = calc_tb_toa(tran_rosen,tbup_rosen,tbdw_rosen,t_surf,emiss,space_tb)
        tb_toa_rosen = np.reshape(tb_toa_rosen,(721,1440))
        tb_maps_rosen[side_band,ifreq,:,:] = tb_toa_rosen

        fig3,ax3 = plot_global_map(np.flipud(tb_toa_rosen-tb_toa_rss),plt_colorbar=True,vmin=-0.5,vmax=0.5,cmap='BrBG',title=f'{freq:.2f} Rosenkranz - RSS')
        fig3.patch.set_facecolor('white')
        fig3.patch.set_alpha(1.0)
        plt.show()
        print

# fig,ax = plt.subplots(2,1)
# ax[0].plot(tb_toa_rss,label='rss 2022')
# ax[0].plot(tb_toa_rosen, label='rosenkranz 2017')
# ax[0].legend()
# ax[1].plot(tb_toa_rss-tb_toa_rosen,label='rosenkranz 2017 - rss 2022')
# ax[1].legend()
# plt.show()
print

