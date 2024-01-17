# README for atm_abs_python_fortran

## Project goals

The final goal of this project is to make an easy-to-use tool to convert model output to MSU/AMSU equivalent top of the atmosphere radiances expressed as antenna temperatures.  The conversion could be performed as several levels of detail.

To make it easy to use, we are developing it in python/xarray.  Ideally it would be able to be dask compatible to enable chunking out the work to multiple processors.

It is important to make the conversion as fast as possible without sacrificing too much accuracy.  To speed things up,
in some cases we wrap previously existing FORTRAN code with python, using the f2py tools that are included in numpy.  f2py really only works for linux environments, so these packages can't built on windows (but maybe one could build a windows wheel using linux)


## Current Status

This project currently contains python wrappers for various atmospheric absorption routines written in Fortran.  The fortran code is included in the project.  The python wrappers are generated using f2py. The build is performed using meson, which probably needs to be installed before the build will work.


## Installation

### Requires:
* numpy
* matplotlib
* xarray
* build

In the top level directory, run
```
python -m build .
```
This will build a wheel for your distribution of python in the /dist directory.

To install the wheel
```
python -m pip install dist/atm_abs-0.0.1-cp310-cp310-linux_x86_64.whl
```
(you will have to modify the second command to match the name of the wheel that was made by the build command if you are using a different version of python.)

## Routines
Most useful routines:

* abs_h2o_rss_2022(t,p,pv,f,h2oabs_rss_out)
* abs_o2_rosen_2017(t,p,pv,f,o2abs_rosen_out)
* abs_o2_rss_2022(t,p,pv,f,o2abs_rss_out)

* atm_tran_multiple_profiles
* fdcldabs(freq,t,rhol,al)
* goff_gratch_rss_2022

The below are not currently used and/or in development:
* fdcldabs_2d_ql
* get_atm_components_rosen_2017
* get_atm_components_rosen_2017_2d
* get_atm_components_rss_2022
* get_atm_components_rss_2022_2d
* trilinear_interpolation


Now make the absorption tables that are used in the interpolation routines.

```
cd make_tables
mkdir tables
python make_absorption_table_cld.AMSU.multiband.py
python make_absorption_table_q.AMSU.multiband.py
```

This confirms that the atm_rtm package as built sucessfully.  The first python command executed in a few seconds.
The second one takes much longer.
If you don't want to wait for this, you can download the tables (for only the RSS Oxygen model) from zenodo:

https://zenodo.org/records/10235106


## Usage

Code in the 'test' directory should show how to use the python wrappers.  The fortran code is in the 'src' directory.

