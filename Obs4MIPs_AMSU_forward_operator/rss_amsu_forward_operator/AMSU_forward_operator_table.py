
from dataclasses import dataclass
from typing import Tuple, Dict
import os
import numpy as np
import calc_tb_sounder.calc_tb_sounder as calc_tb_sounder
from numpy.typing import NDArray
import numpy as np

from .AMSU_constants import AMSU_NOM_EIAS,AMSU_VIEW_ANGLES
from .AMSU_constants import AMSU_A_CHAN5_TLT_WTS,AMSU_A_CHAN5_TMT_WTS,AMSU_A_CHAN7_TTS_WTS,AMSU_A_CHAN9_TLS_WTS



@dataclass
class AMSUForwardOperatorTable:

    AMSU_channel: int
    OxygenAbs_index: int = 4
    VaporAbs_index: int = 4

    def __post_init__(self) -> None:
        self.AMSU_channel = int(self.AMSU_channel)
        self.OxygenAbs_index = int(self.OxygenAbs_index)
        self.VaporAbs_index = int(self.VaporAbs_index)
        self.path_to_data = os.path.join(os.path.dirname(calc_tb_sounder.__file__),'data')
        print(f"Initializing AMSU tables for channel {self.AMSU_channel} from netcdf path: {self.path_to_data}")
        err = calc_tb_sounder.rtm_tables_amsu.read_abs_table_q_amsu_netcdf(self.AMSU_channel,self.path_to_data,self.VaporAbs_index,self.OxygenAbs_index)
        if err != 0:
            raise RuntimeError(f"AMSU table initialization failed (error={err}).")
        print(f"Initializing AMSU cloud_tables for channel {self.AMSU_channel} from netcdf path: {self.path_to_data}")
        err = calc_tb_sounder.rtm_tables_amsu.read_cld_abs_table_amsu_netcdf(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU cloud absorption table initialization failed (error={err}).")
        print(f"Initializing AMSU ocean_emiss_tables for channel {self.AMSU_channel} from netcdf path: {self.path_to_data}")
        err = calc_tb_sounder.rtm_tables_amsu.read_ocean_emiss_table_amsu_netcdf(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU ocean emissivity table initialization failed (error={err}).")
        
        err = calc_tb_sounder.rtm_tables_amsu.read_sea_ice_emiss_table_amsu_netcdf(self.AMSU_channel,self.path_to_data)
        if err != 0:
            raise RuntimeError(f"AMSU sea ice emissivity table initialization failed (error={err}).")

    def calc_vs_fov(
                self,
                t: np.ndarray,
                p: np.ndarray,
                q: np.ndarray,
                cld: np.ndarray,
                cld_mix: int,
                surf_t: np.ndarray,
                surf_p: np.ndarray,
                surf_q: np.ndarray,
                surf_cld: np.ndarray,
                num_levels: int,
                theta: np.ndarray,
                num_views: int,
                wind: np.ndarray,
                land_emiss_in: np.ndarray,
                amsu_channel: int | None = None,
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        """Compute brightness temperatures.

        Notes:
            Inputs t/p/q/cld are expected to be sized ``num_levels + 1`` to
            match the Fortran 0-based vertical indexing.
        """
        need_flip = False
        if p[0] > p[1]:
             need_flip = True
        # The rtm expects the profiles with the smallest pressure (top of atmosphere) first, 
        # but some model data may have the opposite ordering. We check the first two levels 
        # of the input pressure profile to determine if we need to flip the vertical dimension 
        # before passing to the Fortran code. We will flip back the output Tbs to match the 
        # original input ordering.

        amsu_channel_value = self.AMSU_channel if amsu_channel is None else int(amsu_channel)
        num_levels_value = int(num_levels)
        num_views_value = int(num_views)

        t,original_shape_t = self._flatten_3d_to_2d(self._strip_first_dims(t))
        if need_flip:
            t = np.flip(t, axis=0)
        t = np.transpose(t)  # transpose to match Fortran column-major order
        if need_flip:
            p = np.flip(p, axis=0)
        try:
            assert(np.all(p[0:-1]-p[1:]<0.0))
        except AssertionError:
            raise RuntimeError("Pressure profile is not strictly increasing.")
        
        q,_ = self._flatten_3d_to_2d(self._strip_first_dims(q))
        if need_flip:
            q = np.flip(q, axis=0)
        q = np.transpose(q)  # transpose to match Fortran column-major order
        cld,_ = self._flatten_3d_to_2d(self._strip_first_dims(cld))
        if need_flip:
            cld = np.flip(cld, axis=0)
        cld = np.transpose(cld)  # transpose to match Fortran column-major order

        try:
            assert(t.shape[0] == p.shape[0] == q.shape[0] == cld.shape[0] == num_levels_value)
        except AssertionError:
            raise RuntimeError(f"Vertical dimension of t/p/q/cld ({t.shape[0]}/{p.shape[0]}/{q.shape[0]}/{cld.shape[0]}) does not match num_levels ({num_levels_value}).")
        try:
            assert(t.shape[1] == q.shape[1] == cld.shape[1])
        except AssertionError:
            raise RuntimeError(f"Number of profiles in t/q/cld do not match (t:{t.shape[1]}, p:{p.shape[1]}, q:{q.shape[1]}, cld:{cld.shape[1]}).")
        
        t_surf,_ = self._flatten_2d_to_1d(self._strip_first_dims(surf_t))
        p_surf,_ = self._flatten_2d_to_1d(self._strip_first_dims(surf_p))
        q_surf,_ = self._flatten_2d_to_1d(self._strip_first_dims(surf_q))
        cld_surf,_ = self._flatten_2d_to_1d(self._strip_first_dims(surf_cld))
        wind,_ = self._flatten_2d_to_1d(self._strip_first_dims(wind))
        
        t = np.asfortranarray(t.astype(np.float32))
        p = np.asfortranarray(p.astype(np.float32))
        q = np.asfortranarray(q.astype(np.float32))
        cld = np.asfortranarray(cld.astype(np.float32))
        t_surf = np.asfortranarray(t_surf.astype(np.float32))
        p_surf = np.asfortranarray(p_surf.astype(np.float32))
        q_surf = np.asfortranarray(q_surf.astype(np.float32))
        cld_surf = np.asfortranarray(cld_surf.astype(np.float32))
        wind = np.asfortranarray(wind.astype(np.float32))
        land_emiss_in = np.asfortranarray(land_emiss_in.astype(np.float32))

        emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, err = \
            calc_tb_sounder.calc_tb_multiview_amsu.calc_tb_multiview_table_amsu_multi_profiles(
            #num_profiles_value,  # These values are not needed because the Fortran code can infer them from 
            #num_levels_value,    # the array dimensions, and passing them separately will lead to errors 
            #num_views_value,     # because fortran doesnt expect them
            t,
            p,
            q,
            cld,
            int(cld_mix),
            t_surf,
            p_surf,
            q_surf,
            cld_surf,
            amsu_channel_value,
            theta,
            wind,
            land_emiss_in
        )
        
        num_lats = original_shape_t[0] if len(original_shape_t) > 2 else 1
        num_lons = original_shape_t[1] if len(original_shape_t) > 2 else 1

        # reshape back to lat/lon maps
        shp_emiss = emissivity.shape
        num_surfaces = shp_emiss[0]  # should be 3 (ocean, land, sea ice)
        
        emissivity = emissivity.reshape(num_surfaces,num_views_value,num_lats, num_lons)
        surf_wt = surf_wt.reshape(num_surfaces,num_views_value,num_lats, num_lons)
        space_wt = space_wt.reshape(num_surfaces,num_views_value,num_lats, num_lons)
        tb = tb.reshape(num_surfaces,num_views_value,num_lats, num_lons)
        tb_up = tb_up.reshape(num_views_value,num_lats, num_lons)
        tb_dw = tb_dw.reshape(num_views_value,num_lats, num_lons)
        tau = tau.reshape(num_views_value,num_lats, num_lons)

        return emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, int(err)
    
    def _strip_first_dims(self, arr: np.ndarray) -> np.ndarray:
        
        if arr.shape[0] == 1:
            return arr[0]
        else:
            return arr
        
    def _flatten_3d_to_2d(self, arr: np.ndarray) -> np.ndarray:
        original_shape = arr.shape
        if arr.ndim == 3:
            return arr.reshape(-1, arr.shape[2]), original_shape
        else:
            return arr,original_shape
        
    def _flatten_2d_to_1d(self, arr: np.ndarray) -> np.ndarray:
        original_shape = arr.shape
        if arr.ndim == 2:
            return arr.flatten(),original_shape
        else:
            return arr,original_shape
        
    def _dewpoint_to_specific_humidity(self, td_k, p_hpa):
        """
        Converts dew point (K) and pressure (hPa) to specific humidity (kg/kg).
        
        Parameters:
        td_k  : Dew point temperature in Kelvin (float or numpy array)
        p_hpa : Total atmospheric pressure in hPa/mb (float or numpy array)
        
        Returns:
        q     : Specific humidity in kg/kg
        """
        # 1. Calculate Actual Vapor Pressure (e) using Magnus-Tetens (Bolton 1980)
        # Td - 273.15 converts Kelvin to Celsius for the empirical constants
        e = 6.112 * np.exp(17.67 * (td_k - 273.15) / (td_k - 29.65))
    
        # 2. Calculate Specific Humidity (q)
        # Formula: q = (epsilon * e) / (P - (1 - epsilon) * e)
        # where epsilon is the ratio of gas constants (0.622)
        #epsilon = 0.622
        #1-epsilon = 0.378

        q = (0.622 * e) / (p_hpa - (0.378 * e))
        return q

    def compute_tbs(self, model_data: Dict[str, NDArray[np.float32]]) -> Dict[str, NDArray[np.float32]]:
        
        if self.AMSU_channel not in [5, 7, 9]:
            raise ValueError(f"Unsupported AMSU channel: {self.AMSU_channel}. Supported channels are 5, 7, and 9.")
        
        surf_q = self._dewpoint_to_specific_humidity(model_data['surface_dewpoint'], model_data['surface_pressure'])
        
        emissivity, surf_wt, space_wt, tb, tb_up, tb_dw, tau, err = self.calc_vs_fov(
            t=model_data['temperature'],
            p=model_data['levels'],
            q=model_data['specific_humidity'],
            cld=model_data['liquid_content'],
            surf_t = model_data['skin_temperature'],
            surf_p = model_data['surface_pressure'],
            surf_cld = np.zeros_like(surf_q),  # assuming no cloud liquid content at surface
            surf_q = surf_q,
            cld_mix=0,  
            num_levels=len(model_data['levels']),
            theta=np.array(AMSU_NOM_EIAS, dtype=np.float32),  # only nadir
            num_views=len(AMSU_NOM_EIAS),  # only nadir
            wind=model_data['wind_10m'], 
            land_emiss_in=np.full([len(AMSU_NOM_EIAS)], 0.9, dtype=np.float32),  # assuming land emissivity is 0.9
            amsu_channel=self.AMSU_channel
        )

        ocean_frac = 1 - model_data['land_fraction']
        si_frac = model_data['sea_ice_fraction']*ocean_frac
        ocean_frac = ocean_frac - si_frac
        land_frac = model_data['land_fraction']
    
        # for TB, surf_wt and emissivity, we need to combine the ocean, land, and sea ice values using the respective fractions.
        # OCEAN = 0
        # LAND  = 1
        # SEA_ICE = 2

        tb_combined = np.full_like(tau, np.nan, dtype=np.float32)
        emiss_combined = np.full_like(tau, np.nan, dtype=np.float32)
        surf_wt_combined = np.full_like(tau, np.nan, dtype=np.float32)
        num_views = tau.shape[0]
        for iview in range(num_views):
            tb_combined[iview,:,:] = (tb[0,iview,:,:]*ocean_frac + tb[1,iview,:,:]*land_frac + tb[2,iview,:,:]*si_frac) 
            emiss_combined[iview,:,:] = (emissivity[0,iview,:,:]*ocean_frac + emissivity[1,iview,:,:]*land_frac + emissivity[2,iview,:,:]*si_frac)
            surf_wt_combined[iview,:,:] = (surf_wt[0,iview,:,:]*ocean_frac + surf_wt[1,iview,:,:]*land_frac + surf_wt[2,iview,:,:]*si_frac)
         
        output_dict = {}
        num_lats = tb_combined.shape[1]
        num_lons = tb_combined.shape[2]

        if self.AMSU_channel == 5:
            tlt_wts = np.array(AMSU_A_CHAN5_TLT_WTS)
            tmt_wts = np.array(AMSU_A_CHAN5_TMT_WTS)

            TLT_Tbs = np.tensordot(tb_combined, tlt_wts, axes=([0],[0]))
            TMT_Tbs = np.tensordot(tb_combined, tmt_wts, axes=([0],[0]))
            
            output_dict['tbs_TMT'] = TMT_Tbs
            output_dict['tbs_TLT'] = TLT_Tbs

        elif self.AMSU_channel == 7:
            tts_wts = np.array(AMSU_A_CHAN7_TTS_WTS)
            TTS_Tbs = np.tensordot(tb_combined, tts_wts, axes=([0],[0]))
            output_dict['tbs_TTS'] = TTS_Tbs
            
        elif self.AMSU_channel == 9:
            tls_wts = np.array(AMSU_A_CHAN9_TLS_WTS)
            TLS_Tbs = np.tensordot(tb_combined, tls_wts, axes=([0],[0]))
            output_dict['tbs_TLS'] = TLS_Tbs
            
        else:
            raise ValueError(f'Unsupported AMSU_channel: {self.AMSU_channel}')
        print
        return output_dict
    






