"""
This makes the ocn forcing files for an analytical run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g ae0 -t v0 -r backfill -s continuation -d 2020.01.01 -f ocnA0 -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np

from lo_tools import Lfun, zfun, zrfun

import Ofun_nc_xarray
# defaults
verbose = True
if Ldir['testing']:
    from importlib import reload
    reload(Ofun_nc_xarray)
    reload(zrfun)

# this directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# get grid and S info
G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}

# Write files to NetCDF.

tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
#Ofun_nc_xarray.make_clm_file(temp_out_fn_list, NT, out_fn)

# -----------------------------------------

# associate grid types with dimensions
grid_type_dict = {
        'r2dvar': ('eta_rho', 'xi_rho'),
        'u2dvar': ('eta_u', 'xi_u'),
        'v2dvar': ('eta_v', 'xi_v'),
        'r3dvar': ('s_rho', 'eta_rho', 'xi_rho'),
        'u3dvar': ('s_rho', 'eta_u', 'xi_u'),
        'v3dvar': ('s_rho', 'eta_v', 'xi_v')
        }

# # associate variables with dimensions
# dims_dict = {
#         'zeta': ('zeta_time', 'eta_rho', 'xi_rho'),
#         'ubar': ('v2d_time', 'eta_u', 'xi_u'),
#         'vbar': ('v2d_time', 'eta_v', 'xi_v'),
#         'salt': ('salt_time', 's_rho', 'eta_rho', 'xi_rho'),
#         'temp': ('temp_time', 's_rho', 'eta_rho', 'xi_rho'),
#         'u': ('v3d_time', 's_rho', 'eta_u', 'xi_u'),
#         'v': ('v3d_time', 's_rho', 'eta_v', 'xi_v')
#         }
# # assign attributes to variables
# attrs_dict = {
#         'zeta': {'long_name': 'sea surface height climatology', 'units':'meter'},
#         'ubar': {'long_name': 'vertically averaged u-momentum climatology', 'units':'meter second-1'},
#         'vbar': {'long_name': 'vertically averaged v-momentum climatology', 'units':'meter second-1'},
#         'salt': {'long_name': 'salinity climatology', 'units':'g kg-1'},
#         'temp': {'long_name': 'potential temperature climatology', 'units':'Celsius'},
#         'u': {'long_name': 'u-momentum component climatology', 'units':'meter second-1'},
#         'v': {'long_name': 'v-momentum component climatology', 'units':'meter second-1'}
#         }
# write fields to the Dataset
ds = xr.Dataset()

# make the time vector
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=1)
ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])

NT = len(ot_vec)
NZ = S['N']
NR = G['M']
NC = G['L']

# write fields
V = dict()
V['zeta'] = np.zeros((NT, NR, NC))
V['ubar'] = np.zeros((NT, NR, NC-1))
V['vbar'] = np.zeros((NT, NR-1, NC))
V['salt'] = 30 * np.ones((NT, NZ, NR, NC))
V['temp'] = 10 * np.ones((NT, NZ, NR, NC))
V['u'] = np.zeros((NT, NZ, NR, NC-1))
V['v'] = np.zeros((NT, NZ, NR-1, NC))
# create masks
mr2 = np.ones((NT, NR, NC)) * G['mask_rho'].reshape((1, NR, NC))
mr3 = np.ones((NT, NZ, NR, NC)) * G['mask_rho'].reshape((1, 1, NR, NC))
mu2 = np.ones((NT, NR, NC-1)) * G['mask_u'].reshape((1, NR, NC-1))
mu3 = np.ones((NT, NZ, NR, NC-1)) * G['mask_u'].reshape((1, 1, NR, NC-1))
mv2 = np.ones((NT, NR-1, NC)) * G['mask_v'].reshape((1, NR-1, NC))
mv3 = np.ones((NT, NZ, NR-1, NC)) * G['mask_v'].reshape((1, 1, NR-1, NC))
# apply masks
V['zeta'][mr2==0] = np.nan
V['ubar'][mu2==0] = np.nan
V['vbar'][mv2==0] = np.nan
V['salt'][mr3==0] = np.nan
V['temp'][mr3==0] = np.nan
V['u'][mu3==0] = np.nan
V['v'][mv3==0] = np.nan

# combine into a dict for the Dataset
VV = dict()
for vn in V.keys():
    d = zrfun.get_varinfo(vn)
    dims = (d['time_name'],) + grid_type_dict[d['grid_type']]
    VV[vn] = (dims, V[vn])
    
ds = xr.Dataset(VV)
for vn in V.keys():
    d = zrfun.get_varinfo(vn)
    ds[vn].attrs['units'] = d['units']
    ds[vn].attrs['long_name'] = d['long_name']
    
# time coordinates
ds['ocean_time'] = (('ocean_time',), ot_vec)
ds['ocean_time'].attrs['units'] = Lfun.roms_time_units
ds['ocean_time'].attrs['long_name'] = 'ocean time'

# and save to NetCDF
Enc_dict = {vn:enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
# -----------------------------------------
print('- Write clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_ini.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_ini_file(in_fn, out_fn)
print('- Write ini file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc_xarray.make_bry_file(in_fn, out_fn)
print('- Write bry file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    print(ds)
    ds.close()

# check results
nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
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
