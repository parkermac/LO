"""
Code to process the ecology NetCDF bottle data.

Takes # minutes to run.

"""

import pandas as pd
import numpy as np
import xarray as xr
import gsw
import sys
from time import time as Time

from lo_tools import Lfun, zfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology_nc'
otype = 'bottle'
in_dir = Ldir['data'] / 'obs' / source

testing = True

if testing:
    year_list = [2023]
    # pd.set_option('display.max_rows', 500)
    # pd.set_option('display.max_columns', 500)
    # pd.options.display.width = 0 # auto-detect full display width
else:
    year_list = range(1999,2024)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

tt0 = Time()
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    in_fn =  in_dir / 'bottle_2006_2017_fixed.p'
    in_fn = in_dir / ('MarineWaterProfilesAndNutrientsYear' + ys + '.nc')

    ds = xr.open_dataset(in_fn)

    # These have dimension (stations)
    sta_num = ds.station_number.values    # array with values 1 to ## {*}
    sta_lon = ds.Longitude.values
    sta_lat = ds.Latitude.values
    sta_name = [item.decode() for item in ds.Station.values]

    # These have dimension (profiles)
    pro_sta_num = ds.station_index.values # array with values 1 to ## {*}
    pro_num = ds.profile_index.values     # array with values related to profiles {**}
    pro_date_dti = pd.to_datetime(ds.FieldDate.values)

    # These have dimension (obs)
    obs_pro_num = ds.obs_index.values     # array with values related to profiles {**}
    obs_dti = pd.to_datetime(ds.UTCDatetime.values)
    z = -ds.Depth.values

    # # populate a Dataframe "df" with the data from the NetCDF file
    # # ...
            
    # # Renumber cid to be increasing from zero in steps of one.
    # df = obs_functions.renumber_cid(df)
            
    # # (3) retain only selected variables
    # cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
    #     'CT', 'SA', 'DO (uM)',
    #     'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
    #     'TA (uM)', 'DIC (uM)']
    # this_cols = [item for item in cols if item in df.columns]
    # df = df[this_cols]
        
    # print(' - processed %d casts' % ( len(df.cid.unique()) ))

    # df['cruise'] = None
    
    # if (len(df) > 0) and (testing == False):
    #     # Save the data
    #     df.to_pickle(out_fn)
    #     info_df = obs_functions.make_info_df(df)
    #     info_df.to_pickle(info_out_fn)
        
print('Total time = %d sec' % (int(Time()-tt0)))

