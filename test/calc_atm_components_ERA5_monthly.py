
import numpy as np
import os
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import xarray as xr

#from atm_abs import o2abs_19_sub, absoxy_rss_2022_sub
import atm_rtm

import sys
sys.path.append('/mnt/m/job_access/atm_abs_python_fortran/make_tables')
sys.path.append('/mnt/m/job_access/atm_abs_python_fortran/test')
import msu_constants
import amsu_constants

import matplotlib.pyplot as plt
from rss_plotting.global_map import plot_global_map

import time

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

def load_era5_monthly(*,
                      year,
                      month,
                      path='/mnt/n/data/model/ERA5/monthly'):
    
    short_names = {
        'surface_pressure':('PS','sp'),
        '2m_temperature': ('T2m','t2m'),
        'skin_temperature':('TSkin','skt'),
        '2m_dewpoint_temperature': ('TDew','d2m'),
        '10m_wind_speed': ('W10','si10'),
        'sea_ice_cover': ('SeaIce','siconc'),
        'land_sea_mask': ('land_frac','lsm'),
        'total_column_water_vapour':('TPW','prw'),
        'surface_geopotential' : ('GEO','z')
    }
    path = Path(path)

    era5_monthly = {}
    # load the surface data
    path_2d = path / '2D'  / f'{year:04d}'

    vars_2d = ['surface_pressure',
                '2m_temperature',
                '2m_dewpoint_temperature',
                'skin_temperature',
                '10m_wind_speed',
                'sea_ice_cover',
                'land_sea_mask',
                'surface_geopotential']
    
    for var in vars_2d:
        short_name = short_names[var]
        nc_file = path_2d / f'era5_{short_name[0]}_2D_360_181_{year:04d}_{month:02d}.nc'
        ds = xr.open_dataset(nc_file)
        print(ds.keys())
        print(short_name[1])
        era5_monthly[var] = ds[short_name[1]].values

    # load the 3D data
    path_3d = path / '3D'  / f'{year:04d}'

    short_names = {
        'temperature': ('T','t'),
        'specific_humidity': ('Q','q'),
        'specific_cloud_liquid_water_content': ('CLD','clwc'),
        'specific_rain_water_content': ('RAIN','r'),
        'geopotential': ('G','z'),
    }

    vars_3d = ['specific_humidity',
                'temperature',
                'geopotential',
                'specific_cloud_liquid_water_content']


    for var in vars_3d:
        short_name = short_names[var]
        nc_file = path_3d / f'era5_{short_name[0]}_3D_360_181_{year:04d}_{month:02d}.nc'
        ds = xr.open_dataset(nc_file)
        print(ds.keys())
        print(short_name[1])
        era5_monthly[var] = ds[short_name[1]].values
        era5_monthly['pressure'] = ds['level'].values

    return era5_monthly



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

year = 2021
month = 4

era5_monthly = load_era5_monthly(year=year,month=month)

p = era5_monthly['pressure'].astype(np.float32)
p = np.asfortranarray(np.tile(p,(360,181,1)).T)
t = np.asfortranarray(era5_monthly['temperature'].astype(np.float32)[0,:,:,:])
ql = np.asfortranarray(era5_monthly['specific_cloud_liquid_water_content'].astype(np.float32)[0,:,:,:])
q = np.asfortranarray(era5_monthly['specific_humidity'].astype(np.float32)[0,:,:,:])
z = np.asfortranarray(era5_monthly['geopotential'].astype(np.float32)[0,:,:,:])

p_surf = np.asfortranarray(era5_monthly['surface_pressure'].astype(np.float32)[0,:,:])/100.0
t_surf = np.asfortranarray(era5_monthly['skin_temperature'].astype(np.float32)[0,:,:])
dp_surf = np.asfortranarray(era5_monthly['2m_dewpoint_temperature'].astype(np.float32)[0,:,:])
z_surf = np.asfortranarray(era5_monthly['surface_geopotential'].astype(np.float32)[0,:,:])
h_surf =z_surf*re/(9.80665*re - z_surf)

pv_surf = np.asfortranarray(np.zeros((181*360),dtype=np.float32))

atm_rtm.atm_rtm.goff_gratch_rss_2022(np.reshape(dp_surf,-1),
                                    np.ones((181*360),dtype=np.float32)*100.0,
                                    np.reshape(p_surf,-1),
                                    pv_surf)

# fig0,ax0 = plot_global_map(np.flipud(h_surf),plt_colorbar=True,vmin=-200.0,vmax=1100.,cmap='BrBG',title='h')
# fig0,ax0 = plot_global_map(np.flipud(p_surf),plt_colorbar=True,vmin=970.0,vmax=1030.,cmap='BrBG',title='p_surf')

# h =z*re/(9.80665*re - z)
# fig0,ax0 = plot_global_map(np.flipud(h[-1,:,:]),plt_colorbar=True,vmin=-200.0,vmax=1100.,cmap='BrBG',title='z(0)')
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

#convert from geopotential to height
h =z*re/(9.80665*re - z)

#convert from q to pv
pv = p*q/(0.622+0.378*q)

# reverse the order of the profiles so that the surface is at the bottom
h = np.flipud(h)
p = np.flipud(p)
t = np.flipud(t)
q = np.flipud(q)
pv = np.flipud(pv)
ql = np.flipud(ql)

# add the surface values to the bottom of the profiles
h = np.insert(h,0,h_surf,axis=0)
p = np.insert(p,0,p_surf,axis=0)
t = np.insert(t,0,t_surf,axis=0)
q = np.insert(q,0,0.0,axis=0)
pv = np.insert(pv,0,pv_surf,axis=0)
ql = np.insert(ql,0,0.0,axis=0)

rhol = convert_ql_to_rhol(ql,p,q,t)





# find the cloud absorption 

amsu_frequencies = find_frequencies(amsu_channel)

# al = np.zeros_like(p,dtype=np.float32,order='F')
# al_tot = np.zeros_like(p,dtype=np.float32,order='F')
# for freq in amsu_frequencies:
#     print(f'Processing freq {freq:.2f} GHz')
#     atm_rtm.atm_rtm.fdcldabs_2d_ql(freq,p,t,q,ql,al)
#     al_tot = al_tot + al

# cld_abs = al_tot/len(amsu_frequencies)

tht = 0.0
freq = 53.48
shp = h.shape
nlev = shp[0]
nprofile = shp[1]

# tran_arr_test = np.full((nprofile),np.nan,dtype=np.float32,order='F')
# tbdw_arr_test = np.full((nprofile),np.nan,dtype=np.float32,order='F')
# tbup_arr_test = np.full((nprofile),np.nan,dtype=np.float32,order='F')
# atm_rtm.atm_rtm.get_atm_components_rosen_2017_2d(freq, tht, h, p, t, pv, rhol,tran_arr_test,tbup_arr_test,tbdw_arr_test)

# tran_map_test = np.reshape(tran_arr_test,(181,360))
# tbdw_map_test = np.reshape(tbdw_arr_test,(181,360))
# tbup_map_test = np.reshape(tbup_arr_test,(181,360))

# fig1,ax1 = plot_global_map(np.flipud(tran_map_test),plt_colorbar=True,vmin=0.0,vmax=1.0,cmap='BrBG',title='tran_test')
# fig2,ax2 = plot_global_map(np.flipud(tbdw_map_test),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbdw_test')
# fig3,ax3 = plot_global_map(np.flipud(tbup_map_test),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbup_test')

# fig4,ax4 = plot_global_map(np.flipud(np.reshape(pv_surf,(181,360))),plt_colorbar=True,vmin=0.0,vmax=30.,cmap='BrBG',title='pv_surf')

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

abs1_array = abs1_array/1000.0

tht = 0.0
tran_arr = np.full((nprofile),np.nan,dtype=np.float32,order='F')
tbdw_arr = np.full((nprofile),np.nan,dtype=np.float32,order='F')
tbup_arr = np.full((nprofile),np.nan,dtype=np.float32,order='F')

atm_rtm.atm_rtm.atm_tran_multiple_profiles(tht,t,h,abs1_array.T,tran_arr,tbdw_arr,tbup_arr)

tran_arr = np.reshape(tran_arr,(181,360))
tbdw_arr = np.reshape(tbdw_arr,(181,360))
tbup_arr = np.reshape(tbup_arr,(181,360))

fig1,ax1 = plot_global_map(np.flipud(tran_arr),plt_colorbar=True,vmin=0.0,vmax=1.0,cmap='BrBG',title='tran')
fig2,ax2 = plot_global_map(np.flipud(tbdw_arr),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbdw')
fig3,ax3 = plot_global_map(np.flipud(tbup_arr),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbup')

# try removing the lower levels by hand:
h2 = np.zeros_like(h,dtype=np.float32,order='F')
t2 = np.zeros_like(h,dtype=np.float32,order='F')
abs2 = np.zeros_like(h,dtype=np.float32,order='F')

shp = h.shape
nlev = shp[0]
nprofile = shp[1]
num_lev_to_use = np.zeros(nprofile,dtype=np.int32,order='F')

lower_thickness = np.full((nprofile),np.nan,dtype=np.float32,order='F')
for iprofile in range(0,nprofile):
    h1 = h[:,iprofile]
    ok = h1 > h1[0]

    h2x = np.array([h1[0]])
    h2x = np.append(h2x,h1[ok])
    h2[:,iprofile] = h2x[-1]+np.arange(0,nlev)*0.1
    h2[0:len(h2x),iprofile] = h2x
    

    t2x = np.array([t[0,iprofile]])
    t2x = np.append(t2x,t[ok,iprofile])
    t2[:,iprofile] = t2x[-1]
    t2[0:len(t2x),iprofile] = t2x

    abs2x = np.array([abs1_array[iprofile,0]])
    abs2x = np.append(abs2x,abs1_array[iprofile,ok])
    abs2[:,iprofile,] = abs2x[-1]
    abs2[0:len(abs2x),iprofile,] = abs2x

    num_lev_to_use[iprofile] = len(h2x)
    lower_thickness[iprofile] = h2x[1]-h2x[0]

    if iprofile == 2808:
        print(h2[:,iprofile])
        print(t2[:,iprofile])
        print(abs2[:,iprofile])

tht = 0.0
tran_arr2 = np.full((nprofile),np.nan,dtype=np.float32,order='F')
tbdw_arr2 = np.full((nprofile),np.nan,dtype=np.float32,order='F')
tbup_arr2 = np.full((nprofile),np.nan,dtype=np.float32,order='F')

atm_rtm.atm_rtm.atm_tran_multiple_profiles(tht,t2,h2,abs2,tran_arr2,tbdw_arr2,tbup_arr2)

tran_arr2= np.reshape(tran_arr2,(181,360))
tbdw_arr2 = np.reshape(tbdw_arr2,(181,360))
tbup_arr2 = np.reshape(tbup_arr2,(181,360))

fig1,ax1 = plot_global_map(np.flipud(tran_arr2),plt_colorbar=True,vmin=0.0,vmax=1.0,cmap='BrBG',title='tran2')
fig2,ax2 = plot_global_map(np.flipud(tbdw_arr2),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbdw2')
fig3,ax3 = plot_global_map(np.flipud(tbup_arr2),plt_colorbar=True,vmin=100.0,vmax=270.,cmap='BrBG',title='tbup2')


lower_thickness = np.reshape(lower_thickness,(181,360))
fig,ax = plot_global_map(np.flipud(lower_thickness),plt_colorbar=True,vmin=-100.0,vmax=100.,cmap='BrBG',title='lower_thickness')

fig,ax = plot_global_map(np.flipud(np.reshape(num_lev_to_use,(181,360))),plt_colorbar=True,vmin=34.0,vmax=38.,cmap='plasma',title='num_lev_to_use')   
# amsu_channel = 5
# center_freq_list = [53.48,53.71]
# bandwidth_list = [0.17,0.17]
# n_sub_band = 5

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

