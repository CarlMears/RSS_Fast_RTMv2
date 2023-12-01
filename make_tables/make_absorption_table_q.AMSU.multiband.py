#program make_absorption_table
import datetime
import numpy as np
from pathlib import Path
import xarray as xr

# local imports
import atm_rtm
import amsu_constants

def calc_abs_table(*,channel,AMSU_num_freq,O2_model):
    M_W_AIR = 2.8966e-2
    M_W_H2O  = 1.8015324e-2    # kg/mol
    R_GAS    = 8.3145112
    g        = 9.80665
    c_air = R_GAS/M_W_AIR
    c_h2o = R_GAS/M_W_H2O
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


    num_T = 200
    #num_p = 110
    num_q = 150

    T0 = 140.0
    Delta_T = 0.5
    Delta_P = 5.
    Delta_q = 0.001

    # real(4),dimension(0:200,0:110,0:150)  :: abs_table
    # real(4),dimension(0:200,0:110,0:150)  :: abs_table_per_Pa

    t_values = np.arange(140,340.01,Delta_T)
    p_values1 = np.arange(0,30.0,0.1)
    p_values2 = np.arange(30,100,0.5)
    p_values3 = np.arange(100,200,1.0)
    p_values4 = np.arange(200,1101,2.0)
    p_values = np.concatenate((p_values1,p_values2,p_values3,p_values4),axis=0)
    q_values = np.arange(0,num_q+1)*Delta_q

    abs_table = np.zeros((len(q_values),
                          len(p_values),
                          len(t_values)),dtype=np.float32,order='F')

    abs_table_per_Pa = np.zeros((len(q_values),
                                 len(p_values),
                                 len(t_values)),dtype=np.float32,order='F')

    for band_index in range(0,AMSU_num_band):
        for freq_index in range(0,AMSU_num_freq):
            freq = amsu_freq_arr[freq_index,band_index]
            print(f'Processing Channel {channel}: band_index {band_index} freq index {freq_index} freq {freq:.2f} GHz')
            for t_index,t in enumerate(t_values):
                print(f'Processing T {t:.2f} K')
                for p_index,p in enumerate(p_values): 
 
                    pv = p*q_values/(0.622+0.378*q_values)

                    o2abs = np.full_like(pv,np.nan,dtype=np.float32)
                    h2oabs = np.full_like(pv,np.nan,dtype=np.float32)
                    if O2_model == 'RSS_2022':
                        atm_rtm.atm_rtm.abs_o2_rss_2022(np.full_like(pv,t),
                                                        np.full_like(pv,p),
                                                        pv,
                                                        np.full_like(pv,freq),
                                                        o2abs)
                    elif O2_model == 'Rosenkranz_2017':
                        atm_rtm.atm_rtm.abs_o2_rosen_2017(np.full_like(pv,t),
                                                        np.full_like(pv,p),
                                                        pv,
                                                        np.full_like(pv,freq),
                                                        o2abs)
    
                    atm_rtm.atm_rtm.abs_h2o_rss_2022(np.full_like(pv,t),
                                                     np.full_like(pv,p),
                                                     pv,
                                                     np.full_like(pv,freq),
                                                     h2oabs)
                    if p < 0.001:
                        o2abs = np.zeros_like(pv)
                    total_abs = h2oabs + o2abs
                    abs_table[:,p_index,t_index] = abs_table[:,p_index,t_index] + total_abs
                    
#                   ! find density

                    rho_dry = (p-pv)/(c_air*t)
                    rho_vap = pv/(c_h2o*t)
                    
                    # convert to neper/Pa
                    abs_table_per_Pa[:,p_index,t_index] = (abs_table_per_Pa[:,p_index,t_index] + 
                                                            total_abs*0.001/((rho_dry + rho_vap) * g))
                    abs_table_per_Pa[((rho_dry + rho_vap) * g) == 0.0,p_index,t_index] = 0.0

    abs_table = abs_table/(AMSU_num_freq*AMSU_num_band)
    abs_table_per_Pa = abs_table_per_Pa/(AMSU_num_freq*AMSU_num_band)

    date_created = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    abs_table_xr = xr.Dataset(
                                data_vars = dict(
                                    absorptivity = (["specific_humidity","pressure","temperature"],abs_table,
                                    {'long_name':'Absorption coefficient',
                                        'units':'nepers/Km',
                                        'description':'Absorption coefficient for AMSU-A channel 5',
                                        'source_for_O2_abs': O2_model,
                                        'number_of_freqeuncies_sampled_in_band':AMSU_num_freq,
                                        'comments':'Includes both O2 and H2O contributions'})), 
                                coords=dict(
                                        specific_humidity = (["specific_humidity"],q_values,{'units': 'kg/kg'}),
                                        pressure = (["pressure"],p_values,{'units':'hPa'}),
                                        temperature = (["temperature"],t_values,{'units':'K'})
                                        ),
                                attrs={ 'date_created':date_created,
                                        'creator_name':'Carl Mears, Remote Sensing Systems, Santa Rosa, CA',
                                        'creator_email':'mears@remss.com',
                                }
                                )


    nc_file = Path(f'tables/amsu_{channel:02d}_abs_table_q_per_km_{O2_model}.v4.nc')
    abs_table_xr.to_netcdf(nc_file)

    abs_table_per_Pa_xr = xr.Dataset(
                                data_vars = dict(
                                    absorptivity = (["specific_humidity","pressure","temperature"],abs_table_per_Pa,
                                    {'long_name':'Absorption coefficient',
                                        'units':'nepers/Pa',
                                        'description':'Absorption coefficient for AMSU-A channel 5',
                                        'source':'RSS_2022 RTM',
                                        'comments':'Includes both O2 and H2O contributions'})),   
                                coords=dict(
                                        specific_humidity = (["specific_humidity"],q_values,{'units': 'kg/kg'}),
                                        pressure = (["pressure"],p_values,{'units':'hPa'}),
                                        temperature = (["temperature"],t_values,{'units':'K'})
                                        ),
                                attrs={'date_created':date_created,
                                        'creator_name':'Carl Mears, Remote Sensing Systems, Santa Rosa, CA',
                                        'creator_email':'mears@remss.com'}
                                )

    nc_file = Path(f'tables/amsu_{channel:02d}_abs_table_q_per_Pa_{O2_model}.v4.nc')
    abs_table_per_Pa_xr.to_netcdf(nc_file)

    print('Done')

if __name__ == '__main__':
    for channel in np.arange(4,10):
        AMSU_num_freq = 7
        for O2_model in ['RSS_2022','Rosenkranz_2017']:
            calc_abs_table(channel=channel,AMSU_num_freq=AMSU_num_freq,O2_model=O2_model)

# end program