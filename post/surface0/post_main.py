"""
This is the main program for making a SURFACE subset of the daily output.

It creates a single NetCDF file containing only the surface fields
from the history files in a given day.

Testing on mac:

run post_main.py -gtx cas6_v3_lo8b -ro 2 -r backfill -d 2019.07.04 -job surface0

This is fast but something is messed up with the interpolated velocity.  I don't
trust it.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import post_argfun

Ldir = post_argfun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# imports
pth = Path(__file__).absolute().parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import surf_fun

import netCDF4 as nc
from lo_tools import zrfun
import numpy as np
from time import time
import os

start_time = datetime.now()

print(' - Creating surface file(s) for ' + Ldir['date_string'])
f_string = 'f' + Ldir['date_string']

in_dir = Ldir['roms_out'] / Ldir['gtagex'] / f_string
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']

# ======== Create an output file for EDS ==================
fn_list_raw = os.listdir(in_dir)
fn_list = []
for item in fn_list_raw:
    if 'ocean_his' in item and '.nc' in item:
        fn_list.append(in_dir / item)
fn_list.sort()
# shorten the list to be every 2 hours
#fn_list = fn_list[::2]

# Initialize the multi-file input dataset
in_ds = nc.MFDataset(fn_list)

# Initialize
out_fn = out_dir / 'ocean_surface.nc'
print(' - Writing to: ' + str(out_fn))
out_ds = surf_fun.create_ds(out_fn)
surf_fun.add_dims(in_ds, out_ds)
# Add standard 3D-T, rho-grid, variables
vn_list2t = ['Uwind', 'Vwind', 'ocean_time']
vn_list3t = ['salt', 'temp']
slev = -1
surf_fun.add_fields(in_ds, out_ds, vn_list2t, vn_list3t, slev=slev)
# Add custom variables
# - surface velocity on the rho grid
u0 = in_ds['u'][:, slev, :, :].squeeze()
v0 = in_ds['v'][:, slev, :, :].squeeze()
u00 = u0.data
u00[u0.mask] = 0
v00 = v0.data
v00[v0.mask] = 0
u = np.nan * in_ds['salt'][:, slev, :, :].squeeze()
v = u.copy()
u[:, :, 1:-1] = (u00[:, :, 1:] + u00[:, :, :-1])/2
v[:, 1:-1, :] = (v00[:, 1:, :] + v00[:, :-1, :])/2
u[out_ds['salt'][:].mask] = np.nan
v[out_ds['salt'][:].mask] = np.nan
um = np.ma.masked_where(np.isnan(u), u)
vm = np.ma.masked_where(np.isnan(v), v)
#
vv = out_ds.createVariable('u', float, ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
vv.long_name = 'eastward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = um
#
vv = out_ds.createVariable('v', float, ('ocean_time', 'eta_rho', 'xi_rho'), zlib=True)
vv.long_name = 'northward near-surface velocity'
vv.units = 'meter second-1'
vv.time = 'ocean_time'
vv[:] = vm
# Close output Dataset
out_ds.close()
# Close multi-file input Dataset
in_ds.close()

print('\nContents of extracted box file:')
# check on the results
ds = nc.Dataset(out_fn)
for vn in ds.variables:
    print('%s %s max/min = %0.4f/%0.4f' % (vn, str(ds[vn][:].shape),
        np.nanmax(ds[vn][:].data), np.nanmin(ds[vn][:].data)))
ds.close()


# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success'
else:
   result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
post_argfun.finale(Ldir, result_dict)
