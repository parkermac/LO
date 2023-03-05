"""
Code to process DFO bottle data from the NetCDF version.

Performance: 
"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time
import xarray as xr

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()
in_dir = Ldir['data'] / 'obs' / 'dfo1'

# get the datasets
name_list = ['IOS_CTD_Profiles_0.nc','IOS_CTD_Profiles_1.nc']
for name in name_list:
    fn = in_dir / name
    print('\n'+fn.name)
    ds = xr.open_dataset(fn,decode_times=False)

    # initial exploration
    if True:
        for vn in ds.data_vars:
            try:
                units = ds[vn].units
            except:
                units = ''
            try:
                long_name = ds[vn].long_name
            except:
                long_name = ''
            print('%s: %s [%s]' % (vn,long_name,units))
            
        # time range
        dt_ref = datetime(1970,1,1)
        dt0 = dt_ref + timedelta(seconds=float(ds.time[0].values))
        dt1 = dt_ref + timedelta(seconds=float(ds.time[-1].values))
        print(dt0)
        print(dt1)
            
        """
        RESULT:
        IOS_CTD_Profiles_0.nc
        CNDCST01: Sea Water Electrical Conductivity [S/m]
        DOXMZZ01: Oxygen concentration [umol/kg]
        *DOXYZZ01: Oxygen concentration [mL/L]
        PRESPR01: Pressure [decibar]
        PSALST01: Sea Water Practical Salinity [PSS-78]
        PSALST02: Sea Water Practical Salinity [PSS-78]
        SSALST01: Sea Water Salinity [PPT]
        TEMPS601: Sea Water Temperature [degC]
        TEMPS602: Sea Water Temperature [degC]
        TEMPS901: Sea Water Temperature [degC]
        TEMPS902: Sea Water Temperature [degC]
        TEMPST01: Sea Water Temperature [degC]
        agency:  []
        country:  []
        *depth: Depth [m]
        *event_number:  []
        filename:  []
        geographic_area:  []
        instrument_model:  []
        instrument_serial_number:  []
        instrument_type:  []
        *latitude: Latitude [degrees_north]
        *longitude: Longitude [degrees_east]
        mission_id:  []
        platform:  []
        profile: Profile ID []
        project:  []
        scientist:  []
        *sea_water_practical_salinity: Sea Water Practical Salinity [PSS-78]
        sea_water_pressure: Pressure [dbar]
        *sea_water_temperature: Sea Water Temperature [degC]
        *time: Time [seconds since 1970-01-01T00:00:00Z]
        
        Time ranges:
        
        _Profiles_0:
        1965-07-05 21:50:00
        1993-04-08 22:22:11
        
        _Profiles_1:
        1994-02-08 18:57:26
        2021-02-03 16:37:06
        
        """