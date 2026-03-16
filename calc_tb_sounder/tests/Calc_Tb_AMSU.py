#from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple
import os
import numpy as np
import calc_tb_amsu.calc_tb_amsu as calc_tb_amsu



@dataclass
class Calc_Tb_AMSU:
    """Python wrapper around the AMSU multi-view Tb calculator."""

    AMSU_channel: int

    def __post_init__(self) -> None:
        self.AMSU_channel = int(self.AMSU_channel)
        self.path_to_data = os.path.join(os.path.dirname(calc_tb_amsu.__file__),'data')
        err = calc_tb_amsu.rtm_tables_amsu.read_abs_table_q_amsu(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU table initialization failed (error={err}).")
        err = calc_tb_amsu.rtm_tables_amsu.read_cld_abs_table_amsu(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU cloud absorption table initialization failed (error={err}).")
        err = calc_tb_amsu.rtm_tables_amsu.read_ocean_emiss_table_amsu(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU ocean emissivity table initialization failed (error={err}).")
        err = calc_tb_amsu.rtm_tables_amsu.read_sea_ice_emiss_table_amsu(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU sea ice emissivity table initialization failed (error={err}).")
    
    def calc(
        self,
        t: np.ndarray,
        p: np.ndarray,
        q: np.ndarray,
        cld: np.ndarray,
        cld_mix: int,
        num_levels: int,
        theta: np.ndarray,
        num_views: int,
        wind: float,
        land_emiss_in: np.ndarray,
        amsu_channel: int | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        """Compute brightness temperatures.

        Notes:
            Inputs t/p/q/cld are expected to be sized ``num_levels + 1`` to
            match the Fortran 0-based vertical indexing.
        """

        amsu_channel_value = self.AMSU_channel if amsu_channel is None else int(amsu_channel)
        num_levels_value = int(num_levels)
        num_views_value = int(num_views)

        t = self._as_fortran_array(t, np.float32)
        p = self._as_fortran_array(p, np.float32)
        q = self._as_fortran_array(q, np.float32)
        cld = self._as_fortran_array(cld, np.float32)
        theta = self._as_fortran_array(theta, np.float32)
        land_emiss_in = self._as_fortran_array(land_emiss_in, np.float32)

        self._validate_sizes(t, p, q, cld, num_levels_value, "num_levels")
        self._validate_vector(theta, num_views_value, "theta")
        self._validate_vector(land_emiss_in, num_views_value, "land_emiss_in")

        emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, err = calc_tb_amsu.calc_tb_multiview_amsu.calc_tb_multiview_amsu(
            t,
            p,
            q,
            cld,
            int(cld_mix),
            num_levels_value,
            amsu_channel_value,
            theta,
            num_views_value,
            float(wind),
            land_emiss_in,
        )

        return emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, int(err)

    @staticmethod
    def _as_fortran_array(values: np.ndarray, dtype: np.dtype) -> np.ndarray:
        return np.asfortranarray(values, dtype=dtype)

    @staticmethod
    def _validate_sizes(
        t: np.ndarray,
        p: np.ndarray,
        q: np.ndarray,
        cld: np.ndarray,
        num_levels: int,
        label: str,
    ) -> None:
        expected = num_levels + 1
        for name, arr in ("t", t), ("p", p), ("q", q), ("cld", cld):
            if arr.shape[0] != expected:
                raise ValueError(f"{name} must have length {expected} for {label}={num_levels}.")

    @staticmethod
    def _validate_vector(values: np.ndarray, expected: int, label: str) -> None:
        if values.shape[0] != expected:
            raise ValueError(f"{label} must have length {expected}.")


if __name__ == "__main__":
    # Figure where the data file are stored. 
    

    # Initialize the tb calculator for AMSU channel 5.
    tb_calc = Calc_Tb_AMSU(AMSU_channel=5, path_to_data=path_to_data_files)

    

