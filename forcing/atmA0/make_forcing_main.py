"""
This makes the atm forcing files for an analytical run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g ae0 -r backfill -d 2020.01.01 -f atmA0 -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np
from lo_tools import Lfun, zfun, zrfun

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)

# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# get grid and S info, and some sizes
G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
NR = G['M']; NC = G['L']

# Make the time vector.  Here I just have two time points, at the start
# and end of the day, but you could have more, e.g. hourly.  You would still
# want the total time to just be one day.
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=1)
ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])
NT = len(ot_vec)

# Create fields for the state variables.
vn_list = ['Pair','rain','swrad','lwrad_down','Tair','Qair','Uwind','Vwind']

# For now we just fill everything with zeros and nan's
omat = np.zeros((NT, NR, NC))
# mr2 = np.ones((NT, NR, NC)) * G['mask_rho'].reshape((1, NR, NC))
# omat[mr2==0] = np.nan
# NOTE: I when we tried masking for atm fields ROMS did not like it.

for vn in vn_list:
    out_fn = out_dir / (vn + '.nc')
    out_fn.unlink(missing_ok=True)
    ds = xr.Dataset()
    vinfo = zrfun.get_varinfo(vn)
    tname =  vinfo['time_name']
    dims = (tname,) + vinfo['space_dims_tup']
    # You could intervene here by writing something different than omat.
    ds[vn] = (dims, omat.copy())
    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']
    # time coordinate
    ds[tname] = ((tname,), ot_vec)
    ds[tname].attrs['units'] = Lfun.roms_time_units
    ds[tname].attrs['long_name'] = 'ocean time'
    # and save to NetCDF
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()

def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    print(ds)
    ds.close()

# Check results
nc_list = [item + '.nc' for item in vn_list]
if Ldir['testing']:
    # print info about the files to the screen
    for fn in nc_list:
        print_info(out_dir / fn)
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
