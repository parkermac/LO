"""
This is the main program for making the RIVER forcing file, for the
updated ROMS.

NOTE: This is the first use of the new pre/river1 format for the historical and
climatoloy data files. Otherwise based on riv00.

This is the first code to use the new design where forcing goes into a [gridname] folder.

Test on mac in ipython:

run make_forcing_main.py -g cas6 -r backfill -d 2019.07.04 -f riv01 -test True

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta

from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + date_string) / Ldir['frc']

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun

out_fn = out_dir / 'rivers.nc'
out_fn.unlink(missing_ok=True)

# set up the time index for the record
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime(Ldir['date_string'],dsf) - timedelta(days=2.5)
dt1 = datetime.strptime(Ldir['date_string'],dsf) + timedelta(days=4.5)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
yd_ind = pd.Index(dt_ind.dayofyear)

ot_vec = np.array([Lfun.datetime_to_modtime(item) for item in dt_ind])
NT = len(ot_vec)

S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
N = S['N']

grid_fn = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(grid_fn, only_G=True)

# Load a dataframe with info for rivers to get
if Ldir['gridname'] == 'cas6':
    ctag = 'lo_base'
else:
    print('You need to specify a gridname for this ctag.')
    sys.exit()

ri_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag
ri_df_fn = ri_dir / 'river_info.p'
ri_df = pd.read_pickle(ri_df_fn)

# get historical and climatological data files
Ldir['Hflow_fn'] = ri_dir / 'Data_historical' / ('ALL_flow.p')
Ldir['Cflow_fn'] = ri_dir / 'Data_historical' / ('CLIM_flow.p')
Ldir['Ctemp_fn'] = ri_dir / 'Data_historical' / ('CLIM_temp.p')

# get the list of rivers and indices for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')
if Ldir['testing']:
    gri_df = gri_df.loc[['columbia', 'skagit'],:]
NRIV = len(gri_df)

# associate rivers with ones that have temperature climatology data
ri_df = rivfun.get_tc_rn(ri_df)

# get the flow and temperature data for these days
qt_df_dict = rivfun.get_qt(gri_df, ri_df, dt_ind, yd_ind, Ldir, dt1, days)

# Start Dataset
ds = xr.Dataset()

# Add time coordinate
ds['river_time'] = (('river_time',), ot_vec)
ds['river_time'].attrs['units'] = Lfun.roms_time_units
ds['river_time'].attrs['long_name'] = 'river time'

# Add river coordinate
ds['river'] = (('river',), np.arange(1,NRIV+1))
ds['river'].attrs['long_name'] = 'river runoff identification number'

# Add river names
ds['river_name'] = (('river',), list(gri_df.index))
ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NRIV))
ds[vn] = (dims, Vshape)
ds[vn].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        ds[vn] = (('river',), gri_df.idir.to_numpy())
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NRIV)
        ii = 0
        for rn in gri_df.index:
            if gri_df.loc[rn, 'idir'] == 0:
                X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
            elif gri_df.loc[rn, 'idir'] == 1:
                X_vec[ii] = gri_df.loc[rn, 'col_py']
            ii += 1
        ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NRIV)
        ii = 0
        for rn in gri_df.index:
            if gri_df.loc[rn, 'idir'] == 0:
                E_vec[ii] = gri_df.loc[rn, 'row_py']
            elif gri_df.loc[rn, 'idir'] == 1:
                E_vec[ii] = gri_df.loc[rn, 'row_py'] + 1
            ii += 1
        ds[vn] = (('river',), E_vec)
    ds[vn].attrs['long_name'] = vinfo['long_name']
        

# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NRIV))
rr = 0
for rn in gri_df.index:
    qt_df = qt_df_dict[rn]
    flow = qt_df['final'].values
    Q_mat[:,rr] = flow * gri_df.loc[rn, 'isign']
    rr += 1
ds[vn] = (dims, Q_mat)
ds[vn].attrs['long_name'] = vinfo['long_name']
ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature
for vn in ['river_salt', 'river_temp']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')
    if vn == 'river_salt':
        TS_mat = np.zeros((NT, N, NRIV))
    elif vn == 'river_temp':
        TS_mat = np.nan * np.zeros((NT, N, NRIV))
        rr = 0
        for rn in gri_df.index:
            qt_df = qt_df_dict[rn]
            for nn in range(N):
                TS_mat[:, nn, rr] = qt_df['temperature'].values
            rr += 1
    if np.isnan(TS_mat).any():
        print('Error from riv00: nans in river_temp!')
        sys.exit()
    ds[vn] = (dims, TS_mat)
    ds[vn].attrs['long_name'] = vinfo['long_name']
    ds[vn].attrs['units'] = vinfo['units']
    
# Add biology (see the lineup near the end of fennel_var.h)
bvn_list = ['NO3', 'NH4', 'Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo',
        'TIC', 'TAlk', 'LDeC', 'SDeC', 'Oxyg']
for bvn in bvn_list:
    vn = 'river_' + bvn
    vinfo = zrfun.get_varinfo(vn)
    dims = (vinfo['time'],) + ('s_rho', 'river')
    B_mat = np.nan * np.zeros((NT, N, NRIV))
    rr = 0
    for rn in gri_df.index:
        qt_df = qt_df_dict[rn]
        for nn in range(N):
            B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
        rr += 1
    if np.isnan(B_mat).any():
        print('Error from riv00: nans in B_mat for ' + vn)
        sys.exit()
    ds[vn] = (dims, B_mat)
    ds[vn].attrs['long_name'] = vinfo['long_name']
    ds[vn].attrs['units'] = vinfo['units']

# Save to NetCDF
ds.to_netcdf(out_fn)
ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
