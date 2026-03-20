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

The fortran part calculates surface emissivity and tb (and additional intermediate parameters) for all 15 fields of views (FOVs) and 3 types of surfaces (Land, Water and Sea Ice).  (The AMSU-A instrument uses 30 views but the incidence angles are symmetric about nadir.)  In the python part, the results from the various FOVs are combined to make a product targeted to the satellite producted by various groups, e.g TLT and TMT for AMSU channel 5. Using the land fraction and sea ice fraction variables, the results for the difference surfaces are combined.

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

## Usage
```python

    # intialize the AMSU forward operator for choosen channel
    AMSU_channel=5
    OxygenAbs_index = 5  
    #OxygenAbs_index = 4: Rosenkranz 1992 modified: fdabsoxy_1992_modified
    #ioxyOxygenAbs_index = 5: Tretyakov et al 2005 modified: fdabsoxy_tretyakov_modified 
    amsu_op = AMSUForwardOperatorTable(AMSU_channel=5,OxygenAbs_index=OxygenAbs_index)

    #compute the brightness temperatures for the specified month and year
    brightness_temperatures_5 = amsu_op.compute_tbs(model_data)
```
### Input Data

The input data is a dictionary of GCM model output.  It needs to include these variables.

| name | variable | units | dimensions |
| ---- | -------  | ----- | ---------  |
| temperature | temperature profiles | K | [num_lats,num_lons,num_levels] |
| specific_humidity | humidity profiles | kg/kg | [num_lats,num_lons,num_levels] |
| liquid_content  | cloud profiles | kg/kg | [num_lats,num_lons,num_levels] |
| surface_pressure | surface pressure | hPa | [num_lats,num_lons] |
| surface_temperature | surface temperature | K | [num_lats,num_lons] |
| surface_dewpoint  | surface dewpoint | K | [num_lats,num_lons] |
| skin_temperature | temperature of the surface | K | [num_lats,num_lons] |
| wind_10m | 10m wind speed | m/s | [num_lats,num_lons] |
| land_fraction | could be constant | 0.0-1.0 | [num_lats,num_lons] |
| sea_ice_fraction | sea_ice_fraction |0.0-1.0| [num_lats,num_lons] |
| levels | pressure levels | hPa | [num_levels]

All variables should be np.float32
We are assuming that the data are on fixed pressure levels.
(in CMIP, cloud water is not avaliable on pressure levels - we working on that!)

### Output Data
The output data is a dictionary of np.ndarray of the approriate MSU/AMSU products. For example, for AMSU channel 5, the output will have two products, TMT and TLT.
| name | variable | units | dimensions |
| ---- | -------  | ----- | ---------  |
| TLT | tbs_TLT | K | [num_lats,num_lons] |
| TMT | tbs_TMT | K | [num_lats,num_lons] |


## Quick Start

Run the example from examples folder.  Its calculates monthly AMSU maps from ERA5 output:

```bash
cd examples
python test_AMSU_forward_operator_table.py


# or use the debugger to run it....
```
