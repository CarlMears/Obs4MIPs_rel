import importlib.resources  # need this to initialize RSS_surf_emiss

from typing import Dict
from pathlib import Path
from numpy.typing import NDArray
import numpy as np
from importlib import resources


from .AtmWts_method_1 import AtmWt

#class AMSUForwardOperatorTable
class SounderForwardOperatorSimple:

    def __init__(self, sat: str, channel: str) -> None:
        sat = sat.upper()
        allowed_sats = ['MSU','AMSU']        
        if sat not in allowed_sats:
            raise ValueError(f"Invalid satellite: {sat}. Allowed satellites are: {allowed_sats}")
        
        channel = channel.upper()
        allowed_channels = ['TLT','TMT','TTS','TLS']
        if channel not in allowed_channels:
            raise ValueError(f"Invalid channel: {channel}. Allowed channels are: {allowed_channels}")
        self.channel = channel
        AtmWt_dict = {}

        rtm_data_path = resources.files("calc_tb_sounder_simple") / "data" / "wt_tables" 
        
        AtmWt_dict['ocean'] = AtmWt(channel = channel,surface = 'ocean',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)
        AtmWt_dict['land']  = AtmWt(channel = channel,surface = 'land',sat=sat,RTM_Data_Path=rtm_data_path,verbose=True)

        self.AtmWt_dict = AtmWt_dict
        self.sat = sat
        self.channel = channel

    def _strip_first_dims(self, arr: np.ndarray) -> np.ndarray:
        
        if arr.shape[0] == 1:
            return arr[0]
        else:
            return arr

    def compute_tbs(self, 
                    model_data: Dict[str, NDArray[np.float32]],
                    verbose=True,use_skin_temperature=True) -> Dict[str, NDArray[np.float32]]:

        sea_ice_frac = self._strip_first_dims(model_data['sea_ice_fraction'])   
        land_frac = self._strip_first_dims(model_data['land_fraction'])

        sea_ice_frac[np.isnan(sea_ice_frac)] = 0.0
        sea_ice_frac[sea_ice_frac < -1.0] = 0.0
        not_ocean = land_frac + sea_ice_frac
        not_ocean[not_ocean > 1.0] = 1.0

        #assign AtmWt classes.
        AtmWt_MSU_Ocean = self.AtmWt_dict['ocean']
        AtmWt_MSU_Land  = self.AtmWt_dict['land']

        t = model_data['temperature']
        ps = model_data['surface_pressure']
        levels= model_data['levels']

        if use_skin_temperature:
            ts = model_data['skin_temperature']
        else:
            ts = model_data['2m_temperature']

        # calculate the Tbs and  weights.
        if verbose:
            print('Performing Ocean Calculation')
        tbs_ocean,level_wts_ocean,surface_wts_ocean,space_wts_ocean = \
            AtmWt_MSU_Ocean.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)
        if verbose:
            print('Performing Land Calculation')
        tbs_land,level_wts_land, surface_wts_land, space_wts_land   = \
            AtmWt_MSU_Land.AtmLevelWts(temp_profiles=t,ts = ts, ps=ps, levels=levels)

        # combine the land and ocean results together.
        tbs_combined = not_ocean * tbs_land + (1.0 - not_ocean)*tbs_ocean

        return tbs_combined

    
