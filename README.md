# Obs4MIPs Microwave Sounder Forward Operators

## Purpose

The packages here are intended to convert atmospheric profiles and surface information to 
equivalent radiance from microwave sounder channels and/or products.  The most likely use is to convert general circulation climate model output (e.g. CMIP) to MSU or AMSU satellite products for intercomparion and validation.  Our goals for the project for this project are:
* Easy to install and use
* Fast
* More accurate that previous methods used for this purpose
* Transparent about the methods used

## Two Methods
- A simple level weighting method that includes orography and emissivity differences for ocean and land (calc_tb_sounder_simple).  This method is pure python and used numba to speed it up a little
- A fast RTM based on precomputed absorption tables (calc_tb_sounder and Obs4MIPs_forward_operator).  This method used python-wrapped fortran make it faster.  This method is structure so that it can use different absorption models for various atmospheric compeonents.  Currently, it has two different oxygen absorption models implemented.

