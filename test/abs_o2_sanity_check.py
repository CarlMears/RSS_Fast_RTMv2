
import numpy as np
import os
from pathlib import Path
#from atm_abs import o2abs_19_sub, absoxy_rss_2022_sub
import atm_rtm


import matplotlib.pyplot as plt

frequency = np.linspace(50,70,2001,dtype=np.float32)
temperature = np.full_like(frequency,260.0,dtype=np.float32)

pressure_list = [30.0,100.0,300.0,1000.0]
rel_hum_value= 100.0

for pressure_value in pressure_list:
    pressure = np.full_like(frequency,pressure_value,dtype=np.float32)
    rel_num = np.full_like(frequency,rel_hum_value,dtype=np.float32)
    
    pv = np.full_like(frequency,np.nan,dtype=np.float32)
    atm_rtm.atm_rtm.goff_gratch_rss_2022(temperature,pressure,rel_num, pv)
    
    o2abs_rosen = np.full_like(frequency,np.nan,dtype=np.float32)
    o2abs_rss = np.full_like(frequency,np.nan,dtype=np.float32)
    h2oabs_rss = np.full_like(frequency,np.nan,dtype=np.float32)

    atm_rtm.atm_rtm.abs_o2_rss_2022(temperature,pressure,pv,frequency,o2abs_rss)
    atm_rtm.atm_rtm.abs_h2o_rss_2022(temperature,pressure,pv,frequency,h2oabs_rss)
    atm_rtm.atm_rtm.abs_o2_rosen_2017(temperature,pressure,pv,frequency,o2abs_rosen)

    fig,axs = plt.subplots(2,1)
    axs[0].semilogy(frequency,o2abs_rosen,label='Rosenkranz 2017')
    axs[0].semilogy(frequency,o2abs_rss,label='RSS 2022')
    axs[0].semilogy(frequency,h2oabs_rss,label='RSS 2022 H2O')
    axs[1].plot(frequency,o2abs_rosen/o2abs_rss,label='Rosenkranz/RSS')
    axs[1].set_ylim(0.9,1.1)
    axs[0].set_title(f'Atm Absorption, P = {pressure_value} mb, RH = {rel_hum_value}')
    axs[0].legend()
    axs[1].legend()
    axs[1].set_xlabel('Frequency (GHz)')
    axs[0].set_ylabel('Absorption (nepers/km)')

    plot_path = Path('plots/o2_abs_compare')
    os.makedirs(plot_path,exist_ok=True)
    png_file = plot_path / f'o2_abs_compare_{int(pressure_value)}.png'

    plt.savefig(png_file)
    plt.close()

print