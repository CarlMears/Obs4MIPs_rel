# Calc_Tb_AMSU

Python wrapper for the AMSU multi-view brightness temperature calculator in Fortran.

This has only been tested in Linux.

## Usage

```python
import numpy as np
import calc_tb_sounder.calc_tb_sounder as calc_tb_sounder

AMSU_channel = 5
VaporAbs_index = 4
OxygenAbs_index = 4

# Initialize the method by reading in the absorption tables
path_to_data = resources.files("calc_tb_sounder") / "data" 
print(f"Initializing AMSU tables for channel {AMSU_channel} from netcdf path: {path_to_data}")
err = calc_tb_sounder.rtm_tables_amsu.read_abs_table_q_amsu_netcdf(AMSU_channel,path_to_data,VaporAbs_index,OxygenAbs_index)
if err != 0:
    raise RuntimeError(f"AMSU table initialization failed (error={err}).")
print(f"Initializing AMSU cloud_tables for channel {AMSU_channel} from netcdf path: {path_to_data}")
err = calc_tb_sounder.rtm_tables_amsu.read_cld_abs_table_amsu_netcdf(AMSU_channel,path_to_data)
if err != 0:
    raise RuntimeError(f"AMSU cloud absorption table initialization failed (error={err}).")
print(f"Initializing AMSU ocean_emiss_tables for channel {AMSU_channel} from netcdf path: {path_to_data}")
err = calc_tb_sounder.rtm_tables_amsu.read_ocean_emiss_table_amsu_netcdf(AMSU_channel,path_to_data)
if err != 0:
    raise RuntimeError(f"AMSU ocean emissivity table initialization failed (error={err}).")

err = calc_tb_sounder.rtm_tables_amsu.read_sea_ice_emiss_table_amsu_netcdf(AMSU_channel,path_to_data)
if err != 0:
    raise RuntimeError(f"AMSU sea ice emissivity table initialization failed (error={err}).")

emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, err = \
            calc_tb_sounder.calc_tb_multiview_amsu.calc_tb_multiview_table_amsu_multi_profiles(
            t,                      # temperature profiles, K, (num_levels,num_profiles)
            p,                      # pressure levels, hPa, (num_levels)
            q,                      # specific humidity profiles, kg/kg, (num_levels,num_profiles)
            cld,                    # cloud profiles, units determined by cld_mix, (num_levels,num_profiles)
            int(cld_mix),           #! flag for how to intepret cld -- if gt 0 then mixing ratio
            t_surf,                 # surface temperature (num_profiles)
            p_surf,                 # surface pressure (num_profiles)
            q_surf,                 # specific humidity at the surface (kg/kg)
            cld_surf,               # cloud content at the surface, same units at cloud profiles, (num_profiles)
            amsu_channel,           # AMSU channel as an integer
            theta,                  # incidence angle (num_fovs)
            wind,                   # surface wind speed (num profiles) Used over water
            land_emiss_in           # assume land emissivity -- this is a constant ~0.9
        )

outputs:
- emissivity                        # emissivity, (3,num_fovs,num_profiles)  The first index is surface type (OCEAN, LAND, SEAICE)
- surf_wt                           # surface weight (3,num_fovs,num_profiles)
- space_wt                          # space_wt (3,num_fovs,num_profiles)
- tb                                # brightness temperature (3,num_fovs,num_profiles)
- tb_up                             # upwelling brightness temperature from the atmosphere at TOA (3,num_fovs,num_profiles)
- tb_dw                             # downwelling brightness temperature from the atmosphere at the surface (3,num_fovs,num_profiles)
- tau                               # transmissivity from the surface to TOA (num_fovs,num_profiles)


All real values are 4-byte reals


```


### Installation (the easy way)

The easier way is to use precompiled wheels available on github.  They are available here for python 3.10 to 3.14:

https://github.com/CarlMears/Obs4MIPs_rel/releases/tag/version-0.1.0

Find the version that matches you python version and download it. The "manylinux" wheels are for linux that uses GNU C Library (glibc), including Debian, Ubuntu, almalinux and CentOS/RHEL.  The "musllinux" are for versions that use musl libc, mostly Alpine Linux.

After downloading your wheel, in Bash

```
python -m pip install {the name of your wheel}
```

### Installation (the harder way)
Compiling from source.  
#### Requirements:
- linux.  (It might to be possible to install in other unix flavors (probably easier) or in windows (probably harder), but that is out of scope for this README.md)
- meson (https://mesonbuild.com/)
- a fortran compiler, usually GFortran (https://fortran-lang.org/learn/os_setup/install_gfortran/)
- NetCDF-Fortran (https://docs.unidata.ucar.edu/netcdf-fortran/current/index.html)  Installation of these libraries is somewhat complex and you might have to compile from source.  
- numpy >= 1.26 (https://numpy.org/)

After these are installed, in linux:
- clone the repo from github
- cd to the /calc_tb_sounder subdirectory
```
- python -m build .
(this makes a python "wheel" in the /dist subdirectory.  Typically you get a lot of warnings. These can be ignored in the wheel (.whl) is sucessfully generated.)
- python -m pip install dist/{name of your wheel}
```