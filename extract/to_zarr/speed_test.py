"""
Code to test speed of accessing netcdf vs. zarr arrays.

"""

from time import time
import xarray as xr
from lo_tools import Lfun

Ldir = Lfun.Lstart()

local = False

if local:
    # local version
    fn0 = Ldir['roms_out'] / 'cas7_t2_x11b' / 'f2026.05.01' / 'ocean_his_0001.nc'
    fn1 = '/gscratch/macc/parker/tmp/h_01.zarr'
else:
    # s3 version with direct access
    fn0 = 's3://liveocean-pmacc/tmp/ocean_his_0001.nc'
    fn1 = 's3://liveocean-test/cas7_t2_x11b_zarr/f2026.05.01/ocean_his_0001.zarr'
    fn1alt = 'https://s3.kopah.uw.edu/liveocean-test/cas7_t2_x11b_zarr/f2026.05.01/ocean_his_0001.zarr'
    storage_options = {'client_kwargs': {'endpoint_url': 'https://s3.kopah.uw.edu'}, 'anon': True}

tt0 = time()

if local:
    ds = xr.open_dataset(fn0)
else:
    ds = xr.open_dataset(fn0, engine='h5netcdf', storage_options=storage_options)
    
a = ds.salt[0,:,10,10].values
b = a*a
ds.close()
print('netcdf time = %0.3f' % (time()-tt0))

tt0 = time()
if local:
    ds = xr.open_zarr(fn1)
else:
    #ds = xr.open_zarr(fn1, storage_options=storage_options)
    ds = xr.open_zarr(fn1alt)
a = ds.salt[0,:,10,10].values
b = a*a
ds.close()
print('zarr time = %0.3f' % (time()-tt0))