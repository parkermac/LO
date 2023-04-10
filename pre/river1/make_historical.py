"""
Program to gather or update historical records for rivers.

This is meant to append one year at a time to our historical database of river flow
for a given collection (specified by "ctag").

Useage:

run make_historical
- This will append the most recent full year to the current collection,
  or alert you if more years are required.
"""

from lo_tools import Lfun, zfun
from lo_tools import river_functions as rivf

from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import xarray as xr
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-ctag', type=str, default='lo_base')
parser.add_argument('-y0', type=str, default='')
parser.add_argument('-y1', type=str, default='')
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
ctag = args.ctag
y0 = args.y0
y1 = args.y1
testing = args.testing

if len(y0) == 0:
    year0 = datetime.now().year - 1
else:
    year0 = int(y0)
if len(y1) == 0:
    year1 = year0
else:
    year1 = int(y1)

Ldir = Lfun.Lstart()

# directory to work in
ri_dir = Ldir['LOo'] / 'pre' / 'river1' / ctag

# load a dataframe with info for rivers to get
ri_df_fn = ri_dir / 'river_info.p'
ri_df = pd.read_pickle(ri_df_fn)

# Load the as-run rivers from forcing files to fill in gaps in ec records
roms_fn = 'extraction_2017.01.01_2021.12.31.nc' # copied from LiveOcean_output/river
roms_ds = xr.open_dataset(Ldir['LOo'] / 'pre' / 'river1' / ctag / 'Data_roms' / roms_fn)
# NOTE: this will only fill in missing ec data if it covers the right year. Hence
# you may have to update this with a newer extraction and rename roms_fn.

# location for output
out_dir = ri_dir / 'Data_historical'
Lfun.make_dir(out_dir)

# check on what we have already
out_fn = out_dir / 'ALL_flow.p'
if out_fn.is_file():
    have_old = True
    old_df = pd.read_pickle(out_fn)
    old_df.to_pickle(out_dir / 'ALL_flow_prev.p')
else:
    have_old = False
    print('No existing ALL_flow.p. Creating new one.')

if have_old:
    tt = old_df.index
    prev_year0 = tt.year[0]
    prev_year1 = tt.year[-1]
    print('Current range of ALL_flow.p is %d to %d' % (prev_year0, prev_year1))
    if year0 == prev_year1 + 1:
        print('Attempting to append %d-%d to ALL_flow with %d-%d' %
            (year0,year1,prev_year0,prev_year1))
    else:
        print('\nOops! You may need to append more than one year to get to %d' % (year1))
        print('Please specify the range in the arguments.')
        sys.exit()
            
# whether or not we have old data we will get the new data
dt0 = datetime(year0,1,1)
dt1 = datetime(year1,12,31)
days = (dt0, dt1)

# initialize a dict for the river data (one key per river)
qt_dict = dict()

# get USGS river data
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.usgs):
        print('\n'+(' getting usgs ' + rn).center(60,'-'))
        if rn in ['skokomish', 'hamma']:
            rs, qt = rivf.get_usgs_data_custom(rs, days)
        else:
            rs, qt = rivf.get_usgs_data(rs, days)
        if rs.got_data:
            qt_dict[rn] = qt
            print('  got data')
        else:
            print('  ** no data **')
            
# get ec data, a year at a time
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    Qt = pd.Series(dtype='float64') # initialize a Series to concatenate into
    if pd.notnull(rs.ec):
        for year in range(dt0.year, dt1.year + 1):
            
            print((' ' + str(year) + ' getting ec data ' + rn).center(60,'-'))
            
            # First try to get historical year of data.
            rs, qt = rivf.get_ec_data_historical(rs, year)
            if rs.got_data:
                print('  got data from historical')
            elif not rs.got_data:
                # Next try to get it from the current data source.
                this_days = (datetime(year,1,1), datetime(year,12,31))
                rs, qt = rivf.get_ec_data(rs, this_days)
                if rs.got_data:
                    print('  got data from current')
            else:
                # If all else fails get it from a river forcing extraction, in an
                # xarray Dataset called roms_ds.
                try:
                    qti = pd.date_range(start=dt0+timedelta(days=.5), end=dt1+timedelta(days=.5))
                    qt = roms_ds.transport.sel(time=qti, riv=rn)
                    qt = qt.to_series()
                    print('  got data from ROMS forcing')
                    rs['got_data'] = True
                except Exception as e:
                    qt = ''
                    print('  ** no data **')
                    print(e)
                    rs['got_data'] = False
                    
            if rs.got_data:
                this_qt = qt.copy()
                Qt = pd.concat([Qt, this_qt])
            else:
                pass
                
        if not Qt.empty:
            qt_dict[rn] = Qt
        else:
            print('  ** no data **')
            
# clean up and organize into a DataFrame
dt00 = datetime(year0,1,1,12)
dt11 = datetime(year1,12,31,12)

t = pd.date_range(start=dt00, end=dt11)
rn_list = ri_df.index.to_list()

new_df = pd.DataFrame(index=t, columns=rn_list)
for rn in qt_dict.keys():
    r = qt_dict[rn]
    r = r.sort_index() # the ordering of the historical ec data is off
    r = r.reindex(index=t) # this makes missing dates into nans
    r[r<0] = np.nan # there were a few bad points in willapa 2020
    new_df[rn] = r
# organize columns by mean flow
mean_ser = new_df.mean()
mean_ser = mean_ser.sort_values(ascending=False)
new_df = new_df[mean_ser.index] # biggest rivers first

# save results, concatenating with old if appropriate
if have_old:
    combined_df = pd.concat([old_df, new_df])
    combined_df.to_pickle(out_fn)
else:
    if not testing:
        new_df.to_pickle(out_fn)
    else:
        print('Testing: new_df not saved.')
