# Obs4MIPs Microwave Sounder Forward Operators

## Purpose

The packages here are intended to convert atmospheric profiles and surface information to 
equivalent radiance from microwave sounder channels and/or product.  The most likely use is to convert NPW or general circulation climate model output (e.g. CMIP) to MSU or AMSU satellite products for intercomparion and validation.  Our goals for the project for this project are:
* Easy to install and use
* Fast
* More accurate that previous methods used for this purpose
* Transparent about the methods used

## Two different methods

We provide two different methods
1. A simple vertical weighting method.  This method uses the temperature and geopotential at pressure levels and the surface pressure and surface type (land, ice, or water).  It use the total amount of atmosphere at a location (by the surface pressure), so it automatically account for orography.  But it ignore any changes in moisture in the atmosphere (vapor or clouds) and changes for the ocean emissivity due to changes in wind speed.  This method is coded in python, with "numba" acceleration to speed up the loops.


2. A table-based radiative transfer model (RTM).  In this method, we pre-compute absorption terms in the atmosphere (oxygen, vapor, liquid water) for each channel for the Microwave Sounding Unit (MSU) and the Advanced Microwave Sounding Unit (AMSU) and put the results in tables.  A similar method is used to "tablelize" ocean emissivity (depends on incidence angle, wind speed and temperature) and sea ice (depends on incidence angle and temperature). Using interpolated values from the tables as inputs to the RTM, we get a more accurate estimate of radiance without having to the line-by-line calculated billions of times.  This method is coded in combination of python and fortran, using f2py wrappers.

## Layout

- `/sounder_tbs_simple/` - Method 1
- `/calc_tb_sounder/` - the "rtm engine" for Method 2.
- `/Obs4MIPS_AMSU_forward_operator/` - more user friendly interface for Method 2

## Installation
We recommend installing the code in you local python environment.  
### Installing Method 1.
#### Requirements:
- build (https://build.pypa.io/en/stable/installation.html)
- numpy (https://numpy.org/)
- numba (https://numba.pydata.org/)

In linux:

bash


- cd to the /sounder_tbs_simple directory
- python -m build .

(this makes a python "wheel" in the /dist subdirectory)
- python -m pip install dist/{name of your wheel}'

(the name of the wheel depends on your version of python and numpy)

### Installing Method 2

First, install calc_tb_sounder

#### Requirements:
- linux.  (It might to be possible to install in other unix flavors (probably easier) or in windows (probably harder), but that is out of scope for this README.md)
- meson (https://mesonbuild.com/)
- a fortran compiler, usually GFortran (https://fortran-lang.org/learn/os_setup/install_gfortran/)
- NetCDF-Fortran (https://docs.unidata.ucar.edu/netcdf-fortran/current/index.html)  Installation of these libraries is somewhat complex.  If you are conda user, we recommend using conda to install it. (https://anaconda.org/channels/conda-forge/packages/netcdf-fortran/overview)
- numpy >= 1.26 (https://numpy.org/)

After these are installed, in linux:

- cd to the /calc_tb_sounder subdirectory
- python -m build .

(this makes a python "wheel" in the /dist subdirectory.  Typically you get a lot of warnings. These can be ignored in the wheel (.whl) is sucessfully generated.)
- python -m pip install dist/{name of your wheel}'

(the name of the wheel depends on your version of python and numpy)

second





## Quick Start

Run the example from the module folder:

```bash
cd RSS_AMSU_forward_operator
python examples/test_AMSU_forward_table.py
```
