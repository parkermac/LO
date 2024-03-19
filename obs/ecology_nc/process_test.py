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
fn = in_dir / 'test2002.nc'

ds = xr.open_dataset(fn, decode_times=False)

"""
<xarray.Dataset>
Dimensions:                  (obs: 34015, profiles: 306, stations: 41)
Coordinates:
    DepthInterval            (obs) float32 1.0 1.5 2.0 2.5 ... 143.5 144.0 144.5
    FieldDate                (profiles) |S36 b'2002-01-02' ... b'2002-12-05'
    LongitudeDecimalDegrees  (stations) float32 -122.9 -122.4 ... -122.9 -123.5
    LatitudeDecimalDegrees   (stations) float32 47.09 47.29 47.26 ... 47.4 48.13
    Station                  (stations) |S36 b'BUD005' b'CMB003' ... b'PAH003'
Dimensions without coordinates: obs, profiles, stations
Data variables:
    NO3                      (obs) float32 28.59 nan nan nan ... nan nan nan nan
    Salinity                 (obs) float32 27.35 27.54 27.58 ... 32.44 32.44
    Temp                     (obs) float32 8.87 8.878 8.877 ... 8.305 8.304
    UTCDatetime              (obs) |S36 b'2002-01-02T20:41:03Z' ... b'2002-12...
    obs_index                (obs) int16 1249 1249 1249 1249 ... 1554 1554 1554
    row_size                 (profiles) int16 ...
    station_index            (profiles) int16 22 23 24 25 30 ... 30 31 43 44 45
    CastGUID                 (profiles) |S36 b'01D70998-EBBD-4EEC-B803-27A774...
    profile_index            (profiles) int16 1249 1250 1251 ... 1552 1553 1554
Attributes: (12/17)
    Conventions:                   CF-1.11-draft; ACDD 1-3
    title:                         Vertical in situ sensor and discrete water...
    institution:                   Washington State Department of Ecology
    project:                       Long-term Marine Water Quality Monitoring ...
    geospatial_lat_min:            46.4
    geospatial_lat_max:            49.0
    ...                            ...
    source:                        0.5 m bin averaged in situ sensor profilin...
    date_created:                  2024-03-13T20:50:41Z
    history:                       File created using R software and RNetCDF ...
    featureType:                   timeSeriesProfile
    references:                    https://apps.ecology.wa.gov/publications/s...
    comment:                       1) We adjust sensor data for drift for chl...

NOTES:

Each profile has a profile_index (306 values from 1249 to 1554).

Get the N-th item in profile_index:
pind = ds.profile_index[N]

Get the data for this profile:
salt = ds.Salt[ds.obs_index == pind]

Get the station name for this profile:
sind = ds.station_index[N]
sn = ds.Station[sind-1]

"""


# these have dimension (stations) 41
lon = ds.LongitudeDecimalDegrees.values
lat = ds.LatitudeDecimalDegrees.values
sn = ds.Station.values
sn_list = [item.decode() for item in sn]

# These have dimension (profiles) 306
#
sta_ind = ds.station_index.values
# sta_ind: array with values 1 to 52 BUT we only have 41 stations here
#
pind = ds.profile_index.values
# pind: array with values from 1249 to 1554
#
c = ds.CastGUID.values
cast_id = [item.decode() for item in c]
# cast_id: list of weird values like '0A4854FA-63C2-4D50-A84A-83CA063A8A0F'
#
fd = ds.FieldDate.values
fd_dc = [item.decode() for item in fd]
fdi = pd.to_datetime(fd_dc)
# DatetimeIndex with values like: Timestamp('2002-01-02 00:00:00')

# These have dimension (obs) 34015
# converting the times
tu = ds.UTCDatetime.values
tu_dc = [item.decode() for item in tu]
ti = pd.to_datetime(tu_dc) # a DatetimeIndex
temp = ds.Temp.values
salt = ds.Salinity.values
NO3 = ds.NO3.values
depth = ds.DepthInterval.values
oind = ds.obs_index.values
"""
obs_index is the key to pulling a single profile out of the obs arrays.
Each profile has a unique profile_index, and from this you can find (I think)
the time, lon, lat and station name. Then looking for all the entries in
obs_index with this profile index you get the cast
"""

plt.close('all')

for N in range(0,306,30):

    # pulling out a single profile
    pn = pind[N]
    mask = oind == pn
    dd = depth[mask]
    zz = -dd
    tt = temp[mask]
    ss = salt[mask]
    nn = NO3[mask]


    try:
        # get other info about this cast
        sind = sta_ind[N] # index numbering starts at 1
        this_sn = sn_list[sind-1]
        this_lon = lon[sind-1]
        this_lat = lat[sind-1]
        # Print info for each profile, or a subset, and match it to

        print('\nN=%d pn=%d sind=%d sn=%s (lon,lat)=(%0.1f, %0.1f) FieldDate=%s' %
        (N, pn, sind,
        this_sn, this_lon, this_lat,
        fdi[N].strftime('%Y.%m.%d')))
    except IndexError:
        print('\nIndex out of range')

    # plot a single cast with all its info
    pfun.start_plot(figsize=(11,8))
    fig = plt.figure()

    ax = fig.add_subplot(131)
    ax.plot(tt,zz,'-ob')
    ax.set_title('T')
    ax.set_ylim(np.floor(zz.min()),0)

    ax = fig.add_subplot(132)
    ax.plot(ss,zz,'-or')
    ax.set_title('S')
    ax.set_ylim(np.floor(zz.min()),0)

    ax = fig.add_subplot(133)
    ax.plot(nn,zz,'og')
    ax.set_title('NO3')
    ax.set_ylim(np.floor(zz.min()),0)


    plt.show()
    pfun.end_plot()

