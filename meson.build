project('atm_rtm', 'c', 'fortran',
  version : '0.0.1',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

py = import('python').find_installation(pure: false, modules : ['numpy'])
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()


inc_np = include_directories(incdir_numpy, incdir_f2py)

f2py_interface = custom_target('atm_rtmmodule.c',
  input : ['src/atm_rtmpy.f90'],
  output : ['atm_rtmmodule.c', 'atm_rtm-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'atm_rtm', 
             'only:', 'abs_o2_rosen_2017',
                      'abs_o2_rss_2022',
                      'abs_h2o_rss_2022',
                      'goff_gratch_rss_2022', 
                      'get_atm_components_rss_2022',
                      'get_atm_components_rss_2022_2d',
                      'get_atm_components_rosen_2017',
                      'get_atm_components_rosen_2017_2d',
                      'trilinear_interpolation',
                      'fdcldabs_2d',
                      'atm_tran_multiple_profiles',
                      'fdcldabs',
                      'fdcldabs_2d_ql',
                      ':']
)

py.extension_module('atm_rtm',
   'src/RSS_2022/RSS_2022_ABSORPTION.F90',
   'src/rosenkranz_2017/o2abs_19.f90',
   'src/RSS_2022/goff_gratch.f90',
   'src/weighting_function/trig_degrees.f90',
   'src/weighting_function/RSS_2022_RTM_ATMOSPHERE.F90',
   'src/atm_rtmpy.f90',
  f2py_interface,
  incdir_f2py / 'fortranobject.c',
  fortran_args : ['-fno-automatic', '-std=f2008'],
  include_directories: inc_np,
  dependencies : py_dep,
  install : true
)