# README for atm_abs_python_fortran

This project contains python wrappers for various atmospheric absorption routines written in Fortran.  The fortran code is included in the project.  The python wrappers are generated using f2py.

## Installation
In the top level directory, run
```
python -m build .
python -m pip install dist/atm_abs-0.0.1-cp310-cp310-linux_x86_64.whl
```

## Usage
Code in the 'test' directory should show how to use the python wrappers.  The fortran code is in the 'src' directory.

