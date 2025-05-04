"""
Code to test access speed for kopah.

Results, all for 24 files:

mac: 6 sec
apogee: 15 sec
klone: 18 sec (cpu-g2, standard access)
klone: 96 sec (cpu-g2, s3 access)
klone: 44 sec (cpu-g2, https access)

NOTE: I had to add the files to kopah using --acl-public to get it to work.

"""

import pandas as pd
import xarray as xr
import numpy as np
import fsspec
from lo_tools import Lfun
from time import time

Ldir = Lfun.Lstart()

# test of standard access
#fstr = 'f2017.07.04' # test on mac
fstr = 'f2025.05.03'
in_dir = Ldir['parent'] / 'LO_roms' / 'cas7_t0_x4b' / fstr

a = 0
tt0 = time()
cc = 0
for ii in range(2,26):
    hh = ('0000' + str(ii))[-4:]
    fn = in_dir / ('ocean_his_' + hh + '.nc')
    if False:
        # standard access
        ds = xr.open_dataset(fn)
        a += ds.salt.to_numpy()
        ds.close()
    elif False:
        url = 's3://cas7-t0-x4b/' + fstr + '/ocean_his_' + hh + '.nc'
        print(url)
        fs = fsspec.filesystem('s3', anon=True, endpoint_url='https://s3.kopah.uw.edu')
        ds = xr.open_dataset(fs.open(url))
        a += ds.salt.to_numpy()
        ds.close()
    else:
        fs = fsspec.filesystem('https')
        url = 'https://s3.kopah.orci.washington.edu/cas7-t0-x4b/' + fstr + '/ocean_his_' + hh + '.nc'
        print(url)
        ds = xr.open_dataset(fs.open(url))
        a += ds.salt.to_numpy()
        ds.close()
    cc += 1
print('%0.2f sec for %d files' % (time()-tt0, cc))
