
import numpy as np
import os
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import xarray as xr

#from atm_abs import o2abs_19_sub, absoxy_rss_2022_sub
import atm_rtm

import sys
sys.path.append('/mnt/m/job_access/atm_abs_python_fortran/make_tables')
import msu_constants
import amsu_constants

import matplotlib.pyplot as plt
from rss_plotting.global_map import plot_global_map

from amsu_rtm import AMSU_O2_Q_INTERPOLATOR, AMSU_CLD_INTERPOLATOR

re = 6371.0e3

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

channel = 5
AMSU_num_freq = 7
O2_model_rss = 'RSS_2022'
o2_q_interp_rss = AMSU_O2_Q_INTERPOLATOR(channel, 'RSS_2022')
o2_q_interp_rosen = AMSU_O2_Q_INTERPOLATOR(channel, 'Rosenkranz_2017')
cld_interp = AMSU_CLD_INTERPOLATOR(channel)

era5_file = Path('/mnt/a/data/_access_temp/rtm/era5_levels_2022-04-16.nc')
era5_surf_file = Path('/mnt/a/data/_access_temp/rtm/era5_surface_2022-04-16.nc')

print('Starting Reading ERA5')
ds = xr.open_dataset(era5_file)
ds_surf = xr.open_dataset(era5_surf_file)

#use the first time step as an example
time = 0

# set up the arrays for the RTM.  "asfortranarray" is required for the fortran routines.
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
print('Finished Reading ERA5')


print('Preparing ERA5 data from RTM')

# calculate vapor partial pressure at the surface from the dew point.
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

#convert from geopotential to pysical height in meters
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
print('Starting RTM calculations')

print('Calculating absorption coefficients')
# find the cloud absorption 
cld_abs = cld_interp(rhol,t)

# find the clear sky absorption
o2_vap_abs_rss = o2_q_interp_rss(q,p,t)
o2_vap_abs_rosen = o2_q_interp_rosen(q,p,t)

# combine clear sky and cloud absorption
tot_abs_rss = cld_abs + o2_vap_abs_rss
tot_abs_rosen = cld_abs + o2_vap_abs_rosen

# convert to nepers/m

tot_abs_rss = tot_abs_rss/1000.0
tot_abs_rosen = tot_abs_rosen/1000.0

print('Calculating Trans, Tbdw, Tbup')

nom_eias = amsu_constants.AMSU_NOM_EIAS
num_eias = len(nom_eias)

tht = 0.0

# allocate arrays for the results.  Setting order='F' is required for fortran routines.
tran_arr_rss= np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')
tbup_arr_rss = np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')
tbdw_arr_rss = np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')

tran_arr_rosen= np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')
tbup_arr_rosen = np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')
tbdw_arr_rosen = np.full((nprofile,num_eias),np.nan,dtype=np.float32,order='F')

for eia_index,eia in enumerate(nom_eias):
    print(f'Processing EIA {eia:.2f}')
    # this calls the fortran routine analyze the profiles
    atm_rtm.atm_rtm.atm_tran_multiple_profiles(eia,t,h,
                                               tot_abs_rss,
                                               tran_arr_rss[:,eia_index],
                                               tbdw_arr_rss[:,eia_index],
                                               tbup_arr_rss[:,eia_index])

    atm_rtm.atm_rtm.atm_tran_multiple_profiles(eia,t,h,
                                               tot_abs_rosen,
                                               tran_arr_rosen[:,eia_index],
                                               tbdw_arr_rosen[:,eia_index],
                                               tbup_arr_rosen[:,eia_index])

#reshape back to 721 x 1440 x num_eias maps for plotting
tran_map_rss = np.flipud(np.reshape(tran_arr_rss,(721,1440,num_eias)))
tbup_map_rss = np.flipud(np.reshape(tbup_arr_rss,(721,1440,num_eias)))
tbdw_map_rss = np.flipud(np.reshape(tbdw_arr_rss,(721,1440,num_eias)))

tran_map_rosen = np.flipud(np.reshape(tran_arr_rosen,(721,1440,num_eias)))
tbup_map_rosen = np.flipud(np.reshape(tbup_arr_rosen,(721,1440,num_eias)))
tbdw_map_rosen = np.flipud(np.reshape(tbdw_arr_rosen,(721,1440,num_eias)))


print('Finished RTM calculations')

plot_global_map(tbup_map_rss[:,:,0],title='RSS 2022 Tbup Nadir',cmap='plasma',plt_colorbar=True,vmin=0.0,vmax=300.0)
plot_global_map(tbup_map_rss[:,:,14],title='RSS 2022 Tbup Limb',cmap='plasma',plt_colorbar=True,vmin=0.0,vmax=300.0)

plot_global_map(tran_map_rss[:,:,0],title='RSS 2022 Tran Nadir',cmap='plasma',plt_colorbar=True,vmin=0.0,vmax=1.0)
plot_global_map(tran_map_rss[:,:,14],title='RSS 2022 Tbup Limb',cmap='plasma',plt_colorbar=True,vmin=0.0,vmax=1.0)


print

# tb_maps_rss = np.zeros((2,5,721,1440),dtype=np.float32,order='F')
# tb_maps_rosen = np.zeros((2,5,721,1440),dtype=np.float32,order='F')

# for side_band,center_freq in enumerate(center_freq_list):
#     delta_freq = bandwidth_list[side_band]/n_sub_band
#     k_max = np.floor(n_sub_band/2)
#     freq_list = np.arange(-k_max,k_max+1)*delta_freq + center_freq
#     for ifreq,freq in enumerate(freq_list):
#         print(f'Processing freq {freq:.2f} GHz')
#         theta = 0.0
#         emiss = 0.9
#         space_tb = 2.73

#         num_profiles = h.shape[1]

#         tran_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
#         tbup_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
#         tbdw_rss = np.full((num_profiles),np.nan,dtype=np.float32,order='F')

#         tran_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
#         tbup_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')
#         tbdw_rosen = np.full((num_profiles),np.nan,dtype=np.float32,order='F')

#         atm_rtm.atm_rtm.get_atm_components_rss_2022_2d(
#                     freq, 
#                     theta, 
#                     h, 
#                     p, 
#                     t, 
#                     pv, 
#                     clwc,
#                     tran_rss,
#                     tbup_rss,
#                     tbdw_rss)

#         atm_rtm.atm_rtm.get_atm_components_rosen_2017_2d(
#                     freq, 
#                     theta, 
#                     h, 
#                     p, 
#                     t, 
#                     pv, 
#                     clwc,
#                     tran_rosen,
#                     tbup_rosen,
#                     tbdw_rosen)



        # tb_toa_rss = calc_tb_toa(tran_rss,tbup_rss,tbdw_rss,t_surf,emiss,space_tb)
        # tb_toa_rss = np.reshape(tb_toa_rss,(721,1440))
        # tb_maps_rss[side_band,ifreq,:,:] = tb_toa_rss

        # tb_toa_rosen = calc_tb_toa(tran_rosen,tbup_rosen,tbdw_rosen,t_surf,emiss,space_tb)
        # tb_toa_rosen = np.reshape(tb_toa_rosen,(721,1440))
        # tb_maps_rosen[side_band,ifreq,:,:] = tb_toa_rosen

        # fig3,ax3 = plot_global_map(np.flipud(tb_toa_rosen-tb_toa_rss),plt_colorbar=True,vmin=-0.5,vmax=0.5,cmap='BrBG',title=f'{freq:.2f} Rosenkranz - RSS')
        # fig3.patch.set_facecolor('white')
        # fig3.patch.set_alpha(1.0)
        # plt.show()
        # print

# fig,ax = plt.subplots(2,1)
# ax[0].plot(tb_toa_rss,label='rss 2022')
# ax[0].plot(tb_toa_rosen, label='rosenkranz 2017')
# ax[0].legend()
# ax[1].plot(tb_toa_rss-tb_toa_rosen,label='rosenkranz 2017 - rss 2022')
# ax[1].legend()
# plt.show()
print

