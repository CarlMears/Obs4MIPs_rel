# Calc_Tb_AMSU

Python wrapper for the AMSU multi-view brightness temperature calculator.

## Usage


```python
from Calc_Tb_AMSU import Calc_Tb_AMSU

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
