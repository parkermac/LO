"""
Code to test access speed for kopah.

Results:

mac: 6 sec for 24 files
apogee: 15 sec for 24 files
klone: 18 sec for 24 files (cpu-g2, standard access)

"""

import pandas as pd
import xarray as xr
import numpy as np
import fsspec
from lo_tools import Lfun
from time import time

Ldir = Lfun.Lstart()

# test of standard access
#fstr = 'f2017.07.04'
fstr = 'f2025.05.03'
in_dir = Ldir['parent'] / 'LO_roms' / 'cas7_t0_x4b' / fstr

a = 0
tt0 = time()
cc = 0
for ii in range(2,26):
    hh = ('0000' + str(ii))[-4:]
    fn = in_dir / ('ocean_his_' + hh + '.nc')
    if False:
        ds = xr.open_dataset(fn)
        a += ds.salt.to_numpy()
        ds.close()
    else:
        url = 's3://cas7-t0-x4b/' + fstr + '/ocean_his_' + hh + '.nc'
        fs = fsspec.filesystem('s3', anon=True, endpoint_url='https://s3.kopah.uw.edu')
        ds = xr.open_dataset(fs.open(url))
        a += ds.salt.to_numpy()
        ds.close()
        cc += 1
print('%0.2f sec for %d files' % (time()-tt0, cc))
