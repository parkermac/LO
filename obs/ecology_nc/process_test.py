"""
Code to explore the new NetCDF Ecology Archive.
"""

from lo_tools import Lfun
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

in_dir = Ldir['data'] / 'obs' / 'ecology_nc'
# fn = in_dir / '1999to2023CTDnuts.nc'
fn = in_dir / 'test3years.nc'

ds = xr.open_dataset(fn, decode_times=False)

# these have dimension (station) 76
#
lon = ds.LongitudeDecimalDegrees.values
lat = ds.LatitudeDecimalDegrees.values
#
sn = ds.Station.values
sta_name = [item.decode() for item in sn]

# These have dimension (profile) 9350
#
sta_ind = ds.station_index.values
#
c = ds.CastGUID.values
cast_id = [item.decode() for item in c]
#
fd = ds.FieldDate.values
fd_dc = [item.decode() for item in fd]
fdi = pd.to_datetime(fd_dc)

for ii in range(len(sta_ind))[::10]:
# for ii in range(100):
    sind = sta_ind[ii]
    if sind > 0:
        print('%s %s %2d  %s %d' % (cast_id[ii], sta_name[sind - 1],
            sind, str(fdi[sind - 1]), fdi.year[sind-1]))
    else:
        print('ii = %d, sind = %d' % (ii,sind))

# These have dimension (obs) 1218649
#
# converting the times
tu = ds.UTCDatetime.values
tu_dc = [item.decode() for item in tu]
ti = pd.to_datetime(tu_dc)
#
temp = ds.Temp.values
salt = ds.Salinity.values
NO3 = ds.NO3.values
oxygen = ds.DOAdjusted.values
depth = ds.DepthInterval.values

# pulling out a single profile
pn = 10
pind = ds.profile_index.values
mask = pind == pn
dd = depth[mask]
tt = temp[mask]
ss = salt[mask]
nn = NO3[mask]
ox = oxygen[mask]

# df = pd.DataFrame(index=ti, data={'c':c_dc,'t':t,'p':p})
