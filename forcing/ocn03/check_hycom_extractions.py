"""
One-off code to look at what is in our hycom extractions.
"""

from lo_tools import Lfun, zfun
import xarray as xr
import pandas as pd
from lo_tools import hycom_functions as hfun
from datetime import datetime, timedelta
import numpy as np

Ldir = Lfun.Lstart()

Ldir['date_string'] = '2024.12.13'
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

nd_f = np.ceil(Ldir['forecast_days'])
dt0 = this_dt - timedelta(days=2)      
dt1 = this_dt + timedelta(days=int(nd_f) + 2)
print('dt0 = ' + str(dt0))
print('dt1 = ' + str(dt1))

in_dir = Ldir['LOo'] / 'forcing' / 'cas7' / ('f' + Ldir['date_string']) / 'ocn03' / 'Data'

hkeys = ['ssh','u', 'v','t', 's']

# generate ind_dicts
ind_dicts = dict()
for hkey in hkeys:
    in_fn = in_dir / (hkey + '_tyx.nc')
    # specify the sub region of hycom to extract
    aa = hfun.aa
    # convert to hycom format
    north = aa[3]
    south = aa[2]
    west = aa[0] + 360
    east = aa[1] + 360
    # use the results
    ds = xr.open_dataset(in_fn)
    # find selected indices to use with ncks to extract fields
    t = ds.time.values
    tind = pd.DatetimeIndex(t)
    it0 = np.argwhere(tind==dt0)[0][0]
    it1 = np.argwhere(tind==dt1)[0][0]
    x = ds.lon.values
    y = ds.lat.values
    ix0 = zfun.find_nearest_ind(x,west)
    ix1 = zfun.find_nearest_ind(x,east)
    iy0 = zfun.find_nearest_ind(y,south)
    iy1 = zfun.find_nearest_ind(y,north)
    ds.close()
    ind_dict = {'it0':it0, 'it1':it1, 'ix0':ix0, 'ix1':ix1, 'iy0':iy0, 'iy1':iy1}
    ind_dicts[hkey] = ind_dict

# look at the contents of ind_dicts and the hycom extractions
for hkey in ind_dicts.keys():
    print('\n'+ hkey)
    a = ind_dicts[hkey]
    print('it0=%d, it1=%d, ix0=%d, ix1=%d, iy0=%d, iy1=%d' % (a['it0'], a['it1'], a['ix0'], a['ix1'], a['iy0'], a['iy1'], ))
    # Look at the times in each coordintate file
    in_fn = in_dir / (hkey + '_tyx.nc')
    ds = xr.open_dataset(in_fn)
    dt00 = pd.Timestamp(ds.time[a['it0']].values)
    dt11 = pd.Timestamp(ds.time[a['it1']].values)
    ds.close()
    print('ind_dict:         %s to %s' % (dt00.strftime('%Y.%m.%d %H:%M:%S'), dt11.strftime('%Y.%m.%d %H:%M:%S')))
    # compare to the times in the hycom extractions
    in_fn = in_dir / ('h_' + hkey + '.nc')
    ds = xr.open_dataset(in_fn)
    t = ds.time.values
    tind = pd.DatetimeIndex(t)
    print('hycom extraction: %s to %s' % (tind[0].strftime('%Y.%m.%d %H:%M:%S'), tind[-1].strftime('%Y.%m.%d %H:%M:%S')))
    ds.close()




