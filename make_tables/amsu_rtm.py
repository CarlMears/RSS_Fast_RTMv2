
import datetime
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import xarray as xr


import atm_rtm
import msu_constants
import amsu_constants

import matplotlib.pyplot as plt


class AMSU_O2_Q_INTERPOLATOR:

    def __init__(self, amsu_channel, O2_model):
        self.channel = amsu_channel
        self.O2_model = O2_model
        self.nc_file = Path(f'/mnt/m/job_access/atm_abs_python_fortran/make_tables/tables/amsu_{amsu_channel:02d}_abs_table_q_per_km_{O2_model}.v3.nc')
        self.ds = xr.open_dataset(self.nc_file)
        self.abs_table = self.ds['absorptivity'].values
        self.q = self.ds['specific_humidity'].values
        self.p = self.ds['pressure'].values
        self.t = self.ds['temperature'].values
        self.abs_interpolator = RegularGridInterpolator((self.q, self.p, self.t), self.abs_table, bounds_error=False, fill_value=None)

    def __call__(self, q, p, t):
        return self.abs_interpolator((q, p, t))
    
class AMSU_CLD_INTERPOLATOR:

    def __init__(self, amsu_channel):
        self.channel = amsu_channel
        self.nc_file = Path(f'/mnt/m/job_access/atm_abs_python_fortran/make_tables/tables/amsu_{self.channel:02d}_abs_table_cld_per_km.v3.nc')
        self.ds = xr.open_dataset(self.nc_file)
        self.abs_table = self.ds['absorptivity'].values
        self.rhol = self.ds['cloud_density'].values
        self.t = self.ds['temperature'].values
        self.abs_interpolator = RegularGridInterpolator((self.rhol, self.t), self.abs_table, bounds_error=False, fill_value=None)

    def __call__(self, rhol, t):
        return self.abs_interpolator((rhol, t))