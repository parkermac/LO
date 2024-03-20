"""
Code to explore the new NetCDF Ecology Archive

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
Dimensions:         (obs: 34015, profiles: 306, stations: 41)
Coordinates:
    Depth           (obs) float32 ...
    FieldDate       (profiles) |S20 ...
    Longitude       (stations) float32 ...
    Latitude        (stations) float32 ...
    Station         (stations) |S20 ...
Dimensions without coordinates: obs, profiles, stations
Data variables:
    NO3             (obs) float32 ...
    Salinity        (obs) float32 ...
    Temp            (obs) float32 ...
    UTCDatetime     (obs) |S20 ...
    obs_index       (obs) int16 ...
    row_size        (profiles) int16 ...
    station_index   (profiles) int16 ...
    profile_index   (profiles) int16 ...
    station_number  (stations) int16 ...
Attributes: (12/17)
    Conventions:                   CF-1.11-draft; ACDD 1-3
    title:                         Vertical in situ sensor and discrete water...
    institution:                   Washington State Department of Ecology
    project:                       Long-term Marine Water Quality Monitoring ...
    geospatial_lat_min:            46.4
    geospatial_lat_max:            49.0
    ...                            ...
    source:                        0.5 m bin averaged in situ sensor profilin...
    date_created:                  2024-03-20T14:56:04Z
    history:                       File created using R software and RNetCDF ...
    featureType:                   timeSeriesProfile
    references:                    https://apps.ecology.wa.gov/publications/s...
    comment:                       1) We adjust sensor data for drift for chl...

"""

# These have dimension (stations) 41
sta_num = ds.station_number.values    # array with values 1 to 52 {*}
sta_lon = ds.Longitude.values
sta_lat = ds.Latitude.values
sta_name = [item.decode() for item in ds.Station.values]

# These have dimension (profiles) 306
pro_sta_num = ds.station_index.values # array with values 1 to 52 {*}
pro_num = ds.profile_index.values     # array with values from 1249 to 1554 {**}
pro_date = [item.decode() for item in ds.FieldDate.values]
pro_date_dti = pd.to_datetime(pro_date) # DatetimeIndex

# These have dimension (obs) 34015
obs_pro_num = ds.obs_index.values     # array with values from 1249 to 1554 {**}
obs_dt = [item.decode() for item in ds.UTCDatetime.values]
obs_dti = pd.to_datetime(obs_dt) # a DatetimeIndex
temp = ds.Temp.values
salt = ds.Salinity.values
NO3 = ds.NO3.values
z = -ds.Depth.values

plt.close('all')
for N in range(0,len(pro_num)):

    this_pro_num = pro_num[N]
    this_sta_num = pro_sta_num[N]
    this_dstr = pro_date_dti[N].strftime('%Y.%m.%d')

    # get other info about this cast
    sta_ind = np.argwhere(sta_num==this_sta_num)
    sta_ind = sta_ind[0][0]
    this_sta_name = sta_name[sta_ind]
    this_lon = sta_lon[sta_ind]
    this_lat = sta_lat[sta_ind]

    # inspect selected profiles

    NN = 100
    if np.mod(N,NN)==0:

        print('pro_num=%d sta_num=%2d sta_name=%s (lon,lat)=(%0.1f, %0.1f) FieldDate=%s' %
            (this_pro_num, this_sta_num, this_sta_name, this_lon, this_lat, this_dstr))

        mask = obs_pro_num == this_pro_num
        zz = z[mask]
        tt = temp[mask]
        ss = salt[mask]
        nn = NO3[mask]
        ti = obs_dti[mask]
        # print(ti[0].strftime('%Y.%m.%d')) # Result: matches this_dstr

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

