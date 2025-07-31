"""
This makes the ocn forcing files for ROMS, including the banas-fennel bio fields.

It uses GLORYS (glorys, hereafter) fields instead of HYCOM. Since the glorys fields
are daily averages is greatly simplifies the processing.

Testing:

run make_forcing_main.py -g cas7 -r backfill -d 2017.01.01 -f ocnG00 -test True

run make_forcing_main.py -g cas7 -r forecast -d 2025.06.17 -f ocnG00 -test True

NOTE: the main effect of -test True is that it does the interpolation of glorys
fields to the roms grid much faster, by skipping the nearest-neighbor search
step. This is good for testing, but cannot be used for real applications because
it leaves a lot of unfilled cells.

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
import pandas as pd

from lo_tools import Lfun, zfun, zrfun
    
# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# *** automate when to set add_CTD to True ***
if this_dt == datetime(2012,10,7):
    print('WARNING: adding CTD data to extrapolation!!')
    add_CTD = True

# Set time limits for subsequent processing
if Ldir['run_type'] == 'forecast':
    nd_f = np.ceil(Ldir['forecast_days'])
    dt0 = this_dt
    dt1 = this_dt + timedelta(days=int(nd_f))
    dt_list = pd.date_range(dt0,dt1,freq='D')
elif Ldir['run_type'] == 'backfill':
    dt0 = this_dt
    dt1 = this_dt + timedelta(days=1)
if verbose:
    print('dt0 = ' + str(dt0))
    print('dt1 = ' + str(dt1))
# make the list of datetimes (really a DatetimeIndex)
dt_list = pd.date_range(dt0,dt1,freq='D')
# make a list of datestrings to use with file names
dstr_list = [dt.strftime(Lfun.ds_fmt) for dt in dt_list]
# and a dict relating them
dstr_dict = dict(zip(dt_list,dstr_list))

# Add dye to climatology file making use of zrfun.get_varinfo().
tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
ds = xr.Dataset()    
for vn in V.keys():
    # tt00 = time()
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    tname = vinfo['time_name']
    dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
    ds[vn] = (dims, V[vn])
    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']
    # time coordinate
    ds[tname] = ((tname,), ot_vec)
    ds[tname].attrs['units'] = Lfun.roms_time_units
# add other coordinates
for gtag in ['rho','u','v']:
    ds.coords['lon_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lon_'+gtag])
    ds.coords['lat_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lat_'+gtag])
# add depth and z_rho
ds['h'] = (('eta_rho', 'xi_rho'), G['h'])
ds['z_rho_time'] = (('z_rho_time',), ot_vec)
ds['z_rho_time'].attrs['units'] = Lfun.roms_time_units
ds['z_rho'] = (('z_rho_time', 's_rho', 'eta_rho', 'xi_rho'), VV['z_rho'])
# and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
print('- Add dye to clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()


# if Ldir['start_type'] == 'new':
#     # Write initial condition file if needed
#     tt0 = time()
#     in_fn = out_dir / 'ocean_clm.nc'
#     out_fn = out_dir / 'ocean_ini.nc'
#     out_fn.unlink(missing_ok=True)
#     Ofun2_nc.make_ini_file(in_fn, out_fn)
#     print('- Write ini file: %0.2f sec' % (time()-tt0))
#     sys.stdout.flush()

# # Write boundary file
# tt0 = time()
# in_fn = out_dir / 'ocean_clm.nc'
# out_fn = out_dir / 'ocean_bry.nc'
# out_fn.unlink(missing_ok=True)
# Ofun2_nc.make_bry_file(in_fn, out_fn)
# print('- Write bry file: %0.2f sec' % (time()-tt0))
# sys.stdout.flush()

def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    print(ds)
    ds.close()

# # Check results
# if Ldir['start_type'] == 'new':
#     nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
# else:
#     nc_list = ['ocean_clm.nc', 'ocean_bry.nc']
    
# if Ldir['testing']:
#     # open datasets to have a peek manually
#     dsc = xr.open_dataset(out_dir / 'ocean_clm.nc', decode_times=False)
#     if Ldir['start_type'] == 'new':
#         dsi = xr.open_dataset(out_dir / 'ocean_ini.nc', decode_times=False)
#     dsb = xr.open_dataset(out_dir / 'ocean_bry.nc', decode_times=False)
        
# result_dict['result'] = 'SUCCESS'
# for fn in nc_list:
#     if (out_dir / fn).is_file():
#         pass
#     else:
#         result_dict['result'] = 'FAIL'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
