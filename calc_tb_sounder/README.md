# Calc_Tb_AMSU

Python wrapper for the AMSU multi-view brightness temperature calculator.

## Usage

```python
import numpy as np
import calc_tb_sounder.calc_tb_sounder as calc_tb_sounder

calc_tb = Calc_Tb_AMSU(AMSU_channel=5)

emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, err = calc_tb.calc(
    t=t,
    p=p,
    q=q,
    cld=cld,
    cld_mix=cld_mix,
    num_levels=num_levels,
    amsu_channel=None,
    theta=theta,
    num_views=num_views,
    wind=wind,
    land_emiss_in=land_emiss_in,
)
```

Notes:
- Arrays `t`, `p`, `q`, and `cld` must be sized `num_levels + 1` to match the
  Fortran 0-based indexing used in the solver.
- Initialization loads table files from the paths defined in `rtm_tables_AMSU.f90`.
  Ensure those files are available or update the paths accordingly.


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
