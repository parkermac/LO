"""
This code does a spot check looking fore difference between atm0 and atm00.
"""

import xarray as xr
import numpy as np

import atm_fun as afun

from lo_tools import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v0')

d0 = Ldir['LOo'] / 'forcing' / Ldir['gtag']/ 'f2019.07.04' / 'atm0'
d00 = Ldir['LOo'] / 'forcing' / Ldir['gridname']/ 'f2019.07.04' / 'atm00'

ovl = afun.outvar_list
#ovl = [ovl[0]]
for vn in ovl:
    ds0 = xr.open_dataset(d0 / (vn + '.nc'))
    ds00 = xr.open_dataset(d00 / (vn + '.nc'))
    print('\n=== %s ===' % (vn))
    v0 = ds0[vn][-1,:,:].values
    v00 = ds00[vn][-1,:,:].values
    print('0: max = %0.2f, min = %0.2f (N-nan = %d)' % (np.nanmax(v0), np.nanmin(v0), np.sum(np.isnan(v0))))
    print('00: max = %0.2f, min = %0.2f (N-nan = %d)' % (np.nanmax(v00), np.nanmin(v00), np.sum(np.isnan(v00))))
    ds0.close()
    ds00.close()

