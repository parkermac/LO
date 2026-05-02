"""
Code to test speed of accessing netcdf vs. zarr arrays.

"""

from time import time
import xarray as xr
from lo_tools import Lfun

Ldir = Lfun.Lstart()

fn0 = Ldir['roms_out'] / 'cas7_t2_x11b' / 'f2026.05.01' / 'ocean_his_0001.nc'

#fn1 = 'https://s3.kopah.uw.edu/liveocean-pmacc/LO_roms/cas7_t2_x11b_zarr/f2026.05.01/h_01.zarr'

fn1 = '/gscratch/macc/parker/tmp/h_01.zarr'

tt0 = time()
ds = xr.open_dataset(fn0)
a = ds.salt.values
b = a*a
ds.close()
print('netcdf time = %0.3f' % (time()-tt0))

tt0 = time()
ds = xr.open_zarr(fn1)
a = ds.salt.values
b = a*a
ds.close()
print('zarr time = %0.3f' % (time()-tt0))