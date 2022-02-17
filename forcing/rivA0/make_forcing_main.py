"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_forcing_main.py -g ae0 -t v0 -r backfill -s continuation -d 2020.01.01 -f rivA0 -test True

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + date_string) / Ldir['frc']

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)

out_fn = out_dir / 'rivers.nc'
out_fn.unlink(missing_ok=True)

# Make the time vector.  Here I just have two time points, at the start
# and end of the day, but you could have more, e.g. hourly.  You would still
# want the total time to just be one day.
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=1)
ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])
NT = len(ot_vec)

S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
N = S['N']

grid_fn = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(grid_fn, only_G=True)

# get the list of rivers and indices for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')
NRIV = len(gri_df)

# Start Dataset
ds = xr.Dataset()

# Add time coordinate
ds['river_time'] = (('river_time',), ot_vec)
ds['river_time'].attrs['units'] = Lfun.roms_time_units
ds['river_time'].attrs['long_name'] = 'river time'

# Add river coordinate
ds['river'] = (('river',), np.arange(1,NRIV+1))
ds['river'].attrs['long_name'] = 'river runoff identification number'

# Add river names
ds['river_name'] = (('river',), list(gri_df.index))
ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# for Vtransform = 2, even spacing is a good approximation
Vshape = (1/N) * np.ones((N, NRIV))
ds[vn] = (dims, Vshape)
ds[vn].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        ds[vn] = (('river',), gri_df.idir.to_numpy())
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(Nriv)
        ii = 0
        for rn in gri_df.index:
            if gri_df.loc[rn, 'idir'] == 0:
                X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
            elif gri_df.loc[rn, 'idir'] == 1:
                X_vec[ii] = gri_df.loc[rn, 'col_py']
            ii += 1
        
    ds[vn].attrs['long_name'] = vinfo['long_name']
    


    # # Write Eamp
    # vn = 'tide_Eamp'
    # vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    # dims = ('tide_period',) + vinfo['space_dims_tup']
    # ds[vn] = (dims, Eamp)
    # ds[vn].attrs['units'] = vinfo['units']
    # ds[vn].attrs['long_name'] = vinfo['long_name']
    #
    # # Write all other fields
    # for vn in ['tide_Ephase', 'tide_Cangle', 'tide_Cphase', 'tide_Cmax', 'tide_Cmin']:
    #     vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    #     dims = ('tide_period',) + vinfo['space_dims_tup']
    #     ds[vn] = (dims, omat)
    #     ds[vn].attrs['units'] = vinfo['units']
    #     ds[vn].attrs['long_name'] = vinfo['long_name']

# Save to NetCDF
ds.to_netcdf(out_fn)
ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
