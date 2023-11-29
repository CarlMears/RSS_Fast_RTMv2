#program make_absorption_table
import datetime
import numpy as np
from pathlib import Path
import xarray as xr

import atm_rtm
import msu_constants
import amsu_constants



def calc_cld_abs_table(*,channel,AMSU_num_freq):
    M_W_AIR = 2.8966e-2
    M_W_H2O  = 1.8015324e-2    # kg/mol
    R_GAS    = 8.3145112
    g        = 9.80665
    c_air = R_GAS/M_W_AIR
    c_h2o = R_GAS/M_W_H2O
    rd=287.05
    epsilon=0.622
    rv=rd/epsilon
    one_minus_epsilon=1.0-epsilon
    epsilon_ratio=one_minus_epsilon/epsilon
    
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


    print(amsu_freq_arr)



    num_ql = 300
    Delta_T = 0.5
    Delta_rhol = 0.0005

    # real(4),dimension(0:200,0:110,0:150)  :: abs_table
    # real(4),dimension(0:200,0:110,0:150)  :: abs_table_per_Pa

    t_values = np.arange(140,340.01,Delta_T)
    rhol_values = np.arange(0,num_ql+1)*Delta_rhol

    abs_table = np.zeros((len(rhol_values),
                          len(t_values)),
                          dtype=np.float32,
                          order='F')

    for band_index in range(0,AMSU_num_band):
        for freq_index in range(0,AMSU_num_freq):
            freq = amsu_freq_arr[freq_index,band_index]
            print(f'Processing Channel {channel}: band_index {band_index} freq index {freq_index} freq {freq:.2f} GHz')
            for t_index,t in enumerate(t_values):   
                for rhol_index,rhol in enumerate(rhol_values): 
                    abs_cld = np.zeros((1),dtype=np.float32)
                    abs_cld = atm_rtm.atm_rtm.fdcldabs(freq,t,rhol)
                    abs_table[rhol_index,t_index] += abs_cld

    abs_table = abs_table/(AMSU_num_freq*AMSU_num_band)


    date_created = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    abs_table_xr = xr.Dataset(
                                data_vars = dict(
                                    absorptivity = (["cloud_density","temperature"],abs_table,
                                    {'long_name':'Absorption coefficient',
                                        'units':'nepers/Km',
                                        'description':f'Absorption coefficient for AMSU-A channel {channel}',
                                        'number_of_frequencies_sampled_in_band':AMSU_num_freq})), 
                                coords=dict(
                                        cloud_density = (["cloud_density"],rhol_values,{'units': 'kg/m^3'}),
                                        temperature = (["temperature"],t_values,{'units':'K'})
                                        ),
                                attrs={ 'date_created':date_created,
                                        'creator_name':'Carl Mears, Remote Sensing Systems, Santa Rosa, CA',
                                        'creator_email':'mears@remss.com',
                                }
                                )

    nc_file = Path(f'/mnt/m/job_access/atm_abs_python_fortran/make_tables/tables/amsu_{channel:02d}_abs_table_cld_per_km.v3.nc')
    print(f'Writing to: {nc_file}')
    abs_table_xr.to_netcdf(nc_file)

    print('Done')

if __name__ == '__main__':
    for channel in np.arange(4,10):
        AMSU_num_freq = 7
        for O2_model in ['RSS_2022','Rosenkranz_2017']:
            calc_cld_abs_table(channel=channel,AMSU_num_freq=AMSU_num_freq)

# end program