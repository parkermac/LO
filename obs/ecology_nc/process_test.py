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
# fn = in_dir / 'test3years.nc'
fn = in_dir / 'test2002.nc'

ds = xr.open_dataset(fn, decode_times=False)

"""
<xarray.Dataset>
Dimensions:                  (obs: 34015, profiles: 306, stations: 41)
Coordinates:
    DepthInterval            (obs) float32 ...
    FieldDate                (profiles) |S36 b'2002-01-02' ... b'2002-12-05'
    LongitudeDecimalDegrees  (stations) float32 -122.9 -122.4 -122.4 ... nan nan
    LatitudeDecimalDegrees   (stations) float32 47.09 47.29 47.26 ... 47.4 48.13
    Station                  (stations) |S36 b'BUD005' b'CMB003' ... b'PAH003'
Dimensions without coordinates: obs, profiles, stations
Data variables:
    NO3                      (obs) float32 28.59 nan nan nan ... nan nan nan nan
    Salinity                 (obs) float32 27.35 27.54 27.58 ... 32.44 32.44
    Temp                     (obs) float32 8.87 8.878 8.877 ... 8.305 8.304
    UTCDatetime              (obs) |S36 b'2002-01-02T20:41:03Z' ... b'2002-12...
    row_size                 (profiles) int16 ...
    station_index            (profiles) int16 22 23 24 25 30 ... 30 31 43 44 45
    CastGUID                 (profiles) |S36 b'01D70998-EBBD-4EEC-B803-27A774...
    profile_index            (profiles) int16 ...
Attributes: (12/17)
    Conventions:                   CF-1.11-draft; ACDD 1-3
    title:                         Vertical in situ sensor and discrete water...
    institution:                   Washington State Department of Ecology
    project:                       Long-term Marine Water Quality Monitoring ...
    geospatial_lat_min:            46.4
    geospatial_lat_max:            49.0
    ...                            ...
    source:                        0.5 m bin averaged in situ sensor profilin...
    date_created:                  2024-02-24T00:23:40Z
    history:                       File created using R software and RNetCDF ...
    featureType:                   timeSeriesProfile
    references:                    https://apps.ecology.wa.gov/publications/s...
    comment:                       1) We adjust sensor data for drift for chl...
"""

# these have dimension (stations) 41
lon = ds.LongitudeDecimalDegrees.values
# ** last two values of lon are nan
lat = ds.LatitudeDecimalDegrees.values
sn = ds.Station.values
sta_name = [item.decode() for item in sn]

# These have dimension (profiles) 306
sta_ind = ds.station_index.values
# sta_ind: array with values 1 to 52 BUT we only have 41 stations here
c = ds.CastGUID.values
cast_id = [item.decode() for item in c]
# cast_id: list of weird values like '0A4854FA-63C2-4D50-A84A-83CA063A8A0F'
fd = ds.FieldDate.values
fd_dc = [item.decode() for item in fd]
fdi = pd.to_datetime(fd_dc)
# DatetimeIndex with values like: Timestamp('2002-01-02 00:00:00')

# Print info for each profile, or a subset, and match it to
for ii in range(len(cast_id))[::10]:
    sind = sta_ind[ii]
    try:
        print('%s %s %2d  %s %d' % (cast_id[ii], sta_name[sind - 1],
            sind, str(fdi[ii]), fdi.year[ii]))
    except IndexError:
        print('index out of range')

# These have dimension (obs) 34015
# converting the times
tu = ds.UTCDatetime.values
tu_dc = [item.decode() for item in tu]
ti = pd.to_datetime(tu_dc) # a DatetimeIndex
temp = ds.Temp.values
salt = ds.Salinity.values
NO3 = ds.NO3.values
depth = ds.DepthInterval.values

# pulling out a single profile
pind = ds.profile_index.values
pn = pind[0]
mask = pind == pn
dd = depth[mask]
tt = temp[mask]
ss = salt[mask]
nn = NO3[mask]
