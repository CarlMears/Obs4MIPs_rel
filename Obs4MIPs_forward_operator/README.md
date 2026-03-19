# Obs4MIPs Microwave Sounder Forward Operators

## Purpose

The packages here are intended to convert atmospheric profiles and surface information to 
equivalent radiance from microwave sounder channels and/or products.  The most likely use is to convert general circulation climate model output (e.g. CMIP) to MSU or AMSU satellite products for intercomparion and validation.  Our goals for the project for this project are:
* Easy to install and use
* Fast
* More accurate that previous methods used for this purpose
* Transparent about the methods used

## The method in this subdirectory
This is a table-based radiative transfer model (RTM).  In this method, we pre-compute absorption terms in the atmosphere (oxygen, vapor, liquid water) for each channel for the Microwave Sounding Unit (MSU) and the Advanced Microwave Sounding Unit (AMSU) and put the results in tables.  A similar method is used to "tablelize" ocean emissivity (depends on incidence angle, wind speed and temperature) and sea ice (depends on incidence angle and temperature). Using interpolated values from the tables as inputs to the RTM, we get an accurate estimate of radiance without having to the line-by-line calculations billions of times.  This method is coded in combination of python and fortran, using f2py wrappers.

## Internal dependency
The package here uses the package in /calc_tb_sounder which **should be installed first**.

## Installation
We recommend installing the code in you local python environment.  

First, install calc_tb_sounder

#### Requirements:
- linux.  (It might to be possible to install in other unix flavors (probably easier) or in windows (probably harder), but that is out of scope for this README.md)
- python >= 3.10
- calc_tb_sounder (see the calc_tb_sounder directory in this repo)

### Linux installation
in bash
```bash
- clone this repo to your machine
cd Obs4MIPs_forward_operator
python -m build .
#(this makes a python "wheel" in the /dist subdirectory.)
python -m pip install dist/{name of your wheel}

#(the name of the wheel probably will obs4mips_forward_operator-0.1.0-py3-none-any.whl)
```

## Quick Start

Run the example from examples folder:

```bash
cd examples
python test_AMSU_forward_operator_table.py
# or use the debugger to run it....
```
