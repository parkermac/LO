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
# fn = in_dir / 'test2002.nc'
fn = in_dir / 'MarineWaterProfilesAndNutrientsYear2023.nc'

ds = xr.open_dataset(fn)

for vn in ds.data_vars:
    try:
        print('%s: %s' % (vn, ds[vn].units))
    except AttributeError:
        print('%s: None' % (vn))

"""
These results are for test2002.nc.

<xarray.Dataset>
Dimensions:         (obs: 34015, profiles: 306, stations: 41)
Coordinates:
    Depth           (obs) float32 1.0 1.5 2.0 2.5 ... 143.0 143.5 144.0 144.5
    FieldDate       (profiles) datetime64[ns] 2002-01-02 ... 2002-12-05
    Longitude       (stations) float32 -122.9 -122.4 -122.4 ... -122.9 -123.5
    Latitude        (stations) float32 47.09 47.29 47.26 ... 48.08 47.4 48.13
    Station         (stations) |S6 b'BUD005' b'CMB003' ... b'HCB007' b'PAH003'
Dimensions without coordinates: obs, profiles, stations
Data variables:
    NO3             (obs) float32 28.59 nan nan nan nan ... nan nan nan nan nan
    Salinity        (obs) float32 27.35 27.54 27.58 27.6 ... 32.44 32.44 32.44
    Temp            (obs) float32 8.87 8.878 8.877 8.889 ... 8.304 8.305 8.304
    UTCDatetime     (obs) datetime64[ns] 2002-01-02T20:40:32 ... 2002-12-05T2...
    obs_index       (obs) int16 1249 1249 1249 1249 1249 ... 1554 1554 1554 1554
    row_size        (profiles) int16 ...
    station_index   (profiles) int16 22 23 24 25 30 31 33 ... 25 30 31 43 44 45
    profile_index   (profiles) int16 1249 1250 1251 1252 ... 1551 1552 1553 1554
    station_number  (stations) int16 22 23 24 25 30 31 33 ... 21 36 37 39 42 52
Attributes: (12/17)
    Conventions:                   CF-1.11-draft; ACDD 1-3
    title:                         Vertical in situ sensor and discrete water...
    institution:                   Washington State Department of Ecology
    project:                       Long-term Marine Water Quality Monitoring ...
    geospatial_lat_min:            46.4
    geospatial_lat_max:            49.0
    ...                            ...
    source:                        0.5 m bin averaged in situ sensor profilin...
    date_created:                  2024-03-21T23:19:34Z
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
pro_date_dti = pd.to_datetime(ds.FieldDate.values)

# These have dimension (obs) 34015
obs_pro_num = ds.obs_index.values     # array with values from 1249 to 1554 {**}
obs_dti = pd.to_datetime(ds.UTCDatetime.values)
z = -ds.Depth.values

v_dict = dict()
vn_list = ['Temp','Salinity','DOAdjusted','NO3','NH4']
for vn in vn_list:
    v_dict[vn] = ds[vn].values

plt.close('all')
for N in range(0,len(pro_num)):

    # inspect selected profiles

    NN = 100
    if np.mod(N,NN)==0:

        this_pro_num = pro_num[N]
        this_sta_num = pro_sta_num[N]
        this_dstr = pro_date_dti[N].strftime('%Y.%m.%d')

        # get other info about this cast
        sta_ind = np.argwhere(sta_num==this_sta_num)
        sta_ind = sta_ind[0][0]
        this_sta_name = sta_name[sta_ind]
        this_lon = sta_lon[sta_ind]
        this_lat = sta_lat[sta_ind]

        mask = obs_pro_num == this_pro_num
        zz = z[mask]

        ti = obs_dti[mask]
        this_ti0_dstr = ti[0].strftime('%Y.%m.%d')

        pfun.start_plot(figsize=(16,8))
        fig = plt.figure()
        ncol = len(vn_list)
        ii = 1
        for vn in vn_list:
            ax = fig.add_subplot(1,ncol,ii)
            ax.plot(v_dict[vn][mask],zz,'ob')
            ax.set_title(vn)
            ax.set_ylim(np.floor(zz.min()),0)
            if ii == ncol:
                ax.text(.05,.5,'pro_num=%d\nsta_num=%2d\nsta_name=%s\n(lon,lat)=(%0.1f, %0.1f)\nFieldDate=%s\nUTCDatetime=%s' %
                    (this_pro_num, this_sta_num, this_sta_name, this_lon, this_lat, this_dstr, this_ti0_dstr),
                    transform=ax.transAxes)
            ii += 1

        plt.show()
        pfun.end_plot()

