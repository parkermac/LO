"""
This code does a spot check looking fore difference between riv0 and riv00.
"""

import xarray as xr
import numpy as np

from lo_tools import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v0')

d0 = Ldir['LOo'] / 'forcing' / Ldir['gtag']/ 'f2019.07.04' / 'riv1'
d00 = Ldir['LOo'] / 'forcing' / Ldir['gtag']/ 'f2019.07.04' / 'riv00'

ds0 = xr.open_dataset(d0 / 'rivers.nc')
ds00 = xr.open_dataset(d00 / 'rivers.nc')


for vn in ds0.data_vars:
    if vn in ds00.data_vars:
        print('\n=== %s ===' % (vn))
        v0 = ds0[vn].values
        v00 = ds00[vn].values
        try:
            print('0 : max = %0.2f, min = %0.2f' % (np.nanmax(v0), np.nanmin(v0)))
            print('00: max = %0.2f, min = %0.2f' % (np.nanmax(v00), np.nanmin(v00)))
        except TypeError:
            print('skipping %s' % (vn))

