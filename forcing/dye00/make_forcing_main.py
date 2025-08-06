"""
This adds dye to ocean and river forcing files.

It relies on -gtx to find which forcing files it should add to,
looking in forcing_list.csv.

Testing:

run make_forcing_main.py -g cas7 -r backfill -d 2025.06.17 -f dye00 -gtx cas7_dye0_x11b -ro 1

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

from lo_tools import Lfun, zfun, zrfun, Ofun2_nc
    
# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / Ldir['frc']
# NOTE: we will also make new directories for the forcing we add dye to, by appending
# Ldir['frc'] to the original forcing firectory, so for example ocnG00 would give rise to ocnG00dye00.

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# Look to see if there is a user instance of this dot_in, and if so, use it.
user_dot_in_dir = Ldir['LOu'] / 'dot_in' / Ldir['gtagex']
dot_in_dir = Ldir['LO'] / 'dot_in' / Ldir['gtagex']
if user_dot_in_dir.is_dir():
    dot_in_dir = user_dot_in_dir

# find which forcing files to work on, based on -gtx
force_dict = dict()
with open(dot_in_dir / 'forcing_list.csv', 'r') as f:
    for line in f:
        which_force, force_choice = line.strip().split(',')
        force_dict[which_force] = force_choice

# Add dye to climatology file making use of zrfun.get_varinfo().
tt0 = time()
in_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + Ldir['date_string']) / force_dict['ocn']
in_fn = in_dir /'ocean_clm.nc'
out_dir = in_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
    ('f' + Ldir['date_string']) / (force_dict['ocn']+Ldir['frc'])
Lfun.make_dir(out_dir)
out_fn = out_dir /'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
ds = xr.open_dataset(in_fn)
s = ds.salt # DataArray
vinfo = zrfun.get_varinfo('dye_', vartype='climatology')
for dn in range(2):
    print(dn)
    dn_string = ('000' + str(dn+1))[-2:]
    dye_name = 'dye_' + dn_string
    d = 0 * s
    d.name=dye_name
    ds[dye_name] = d

# for vn in V.keys():
#     # tt00 = time()
#     vinfo = zrfun.get_varinfo(vn, vartype='climatology')
#     tname = vinfo['time_name']
#     dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
#     ds[vn] = (dims, V[vn])
#     ds[vn].attrs['units'] = vinfo['units']
#     ds[vn].attrs['long_name'] = vinfo['long_name']
#     # time coordinate
#     ds[tname] = ((tname,), ot_vec)
#     ds[tname].attrs['units'] = Lfun.roms_time_units
# # add other coordinates
# for gtag in ['rho','u','v']:
#     ds.coords['lon_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lon_'+gtag])
#     ds.coords['lat_'+gtag] = (('eta_'+gtag, 'xi_'+gtag), G['lat_'+gtag])
# # add depth and z_rho
# ds['h'] = (('eta_rho', 'xi_rho'), G['h'])
# ds['z_rho_time'] = (('z_rho_time',), ot_vec)
# ds['z_rho_time'].attrs['units'] = Lfun.roms_time_units
# ds['z_rho'] = (('z_rho_time', 's_rho', 'eta_rho', 'xi_rho'), VV['z_rho'])
# and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
# ds.close()
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

# result_dict['end_dt'] = datetime.now()
# ffun.finale(Ldir, result_dict)
