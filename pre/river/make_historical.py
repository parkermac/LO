"""
Program to gather historical records for rivers.
"""

from lo_tools import Lfun, zfun
from lo_tools import river_functions as rivf

from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import xarray as xr

Ldir = Lfun.Lstart()

# set the gtag for this case (could make it a command line argument)
gtag = 'cas6_v3'
# Load a dataframe with info for rivers to get
ri_fn = Ldir['LOo'] / 'pre' / 'river' / gtag / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

use_roms = True
testing = False
ec_gap_year = 2021 # the year with missing ec data

# Use these to concatenate the results with previous ones
add_to_previous = True
prev_year0 = 1980
prev_year1 = 2020

if use_roms:
    
    # Load the as-run rivers from forcing files to fill in gaps in ec records
    roms_fn = 'extraction_2017.01.01_2021.12.31.nc' # copied from LiveOcean_output/river
    roms_ds = xr.open_dataset(Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_roms' / roms_fn)
    # only needed if we expect a gap between ec_historical and ec current

# location for output
out_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
Lfun.make_dir(out_dir)

# set time range
year0 = 2021
year1 = 2021

if testing:
    ri_df = ri_df.loc[['skagit','fraser']]
        
dt0 = datetime(year0,1,1)
dt1 = datetime(year1,12,31)
days = (dt0, dt1)

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
        if rs['got_data']:
            qt_dict[rn] = qt
        else:
            print(rs)

# get ec data, a year at a time
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    Qt = pd.Series(dtype='float64') # initialize a Series to concatenate into
    if pd.notnull(rs.ec):
        for year in range(dt0.year, dt1.year + 1):
            this_days = (datetime(year,1,1), datetime(year,12,31))
            if year > ec_gap_year:
                print((' ' + str(year) + ' getting ec ' + rn).center(60,'-'))
                rs, qt = rivf.get_ec_data(rs, this_days)
                
            # This section can be used as needed to fill gaps.
            elif use_roms and (year == ec_gap_year):
                try:
                    qti = pd.date_range(start=dt0+timedelta(days=.5), end=dt1+timedelta(days=.5))
                    qt = roms_ds.transport.sel(time=qti, riv=rn)
                    qt = qt.to_series()
                    print((' ' + str(year) + ' getting ec from ROMS forcing ' + rn).center(60,'-'))
                    rs['got_data'] = True
                except Exception as e:
                    qt = ''
                    print(e)
                    rs['got_data'] = False
            
            elif year < ec_gap_year:
                print((' ' + str(year) + ' getting ec historical ' + rn).center(60,'-'))
                rs, qt = rivf.get_ec_data_historical(rs, year)
                
            if rs.got_data:
                this_qt = qt.copy()
                Qt = pd.concat([Qt, this_qt])
            else:
                pass
                
        if not Qt.empty:
            qt_dict[rn] = Qt
            
# clean up and organize into a DataFrame
dt00 = datetime(year0,1,1,12)
dt11 = datetime(year1,12,31,12)

t = pd.date_range(start=dt00, end=dt11)
rn_list = ri_df.index.to_list()

all_df = pd.DataFrame(index=t, columns=rn_list)
for rn in rn_list:
    r = qt_dict[rn]
    r = r.sort_index() # the ordering of the historical ec data is off
    r = r.reindex(index=t) # this makes missing dates into nans
    r[r<0] = np.nan # there were a few bad points in willapa 2020
    all_df[rn] = r
# organize columns by mean flow
mean_ser = all_df.mean()
mean_ser = mean_ser.sort_values(ascending=False)
all_df = all_df[mean_ser.index] # biggest rivers first
new_fn = out_dir / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p')
all_df.to_pickle(out_dir / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p'))

# concatenate with previous historical data
if add_to_previous:
    prev_fn = out_dir / ('ALL_flow_' + str(prev_year0) + '_' + str(prev_year1) + '.p')
    combined_fn = out_dir / ('ALL_flow_' + str(prev_year0) + '_' + str(year1) + '.p')
    prev_df = pd.read_pickle(prev_fn)
    combined_df = pd.concat([prev_df, all_df])
    combined_df.to_pickle(combined_fn)

# PLOTTING
if testing:
    import matplotlib.pyplot as plt
    plt.close('all')

    NP = len(rn_list)
    NR, NC = zfun.get_rc(NP)
    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
    ii = 0
    for rn in qt_dict.keys():
        ir, ic = zfun.get_irc(ii, NC)
        ax = axes[ir, ic]
        all_df[rn].plot(ax=ax)
        ax.set_xlim(dt00, dt11)
        ax.text(.05, .9, rn, transform=ax.transAxes)
        ii += 1

    plt.show()


