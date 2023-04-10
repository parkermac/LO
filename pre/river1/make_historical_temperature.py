"""
Program to gather historical temperature records for rivers.

Patterned off of make_historical, but simpler because we only query
the last requested year of EC data.
"""

from lo_tools import Lfun, zfun
from lo_tools import river_functions as rivf

from datetime import datetime, timedelta
import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-ctag', type=str, default='lo_base')
parser.add_argument('-year0', type=int, default=1980)
parser.add_argument('-year1', type=int, default=2020)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
ctag = args.ctag
year0 = args.year0
year1 = args.year1
testing = args.testing

Ldir = Lfun.Lstart()

# directory to work in
ri_dir0 = Ldir['LOo'] / 'pre' / 'river1' / ctag

# load a dataframe with info for rivers to get
ri_df_fn = ri_dir0 / 'river_info.p'
ri_df = pd.read_pickle(ri_df_fn)

# location for output
out_dir = ri_dir0 / 'Data_historical'
Lfun.make_dir(out_dir)

# check on what we have already
out_fn = out_dir / 'ALL_temperature.p'
if out_fn.is_file():
    have_old = True
    old_df = pd.read_pickle(out_fn)
    old_df.to_pickle(out_dir / 'ALL_temperature_prev.p')
else:
    have_old = False
    print('No existing ALL_temperature.p. Creating new one.')

if testing:
    ri_df = ri_df.loc[['skagit','fraser']]
    year0 = datetime.now().year - 1
    year1 = year0
        
dt0 = datetime(year0,1,1)
dt1 = datetime.now().year - 1 # always go to last full year to get EC data
days = (dt0, dt1)

qt_dict = dict()
# get USGS river data
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.usgs):
        print('\n'+(' getting usgs ' + rn).center(60,'-'))
        rs, qt = rivf.get_usgs_data(rs, days, temperature=True)
        if rs.got_data:
            qt_dict[rn] = qt
        else:
            print(' - no data -')

# get ec data, just for most recent year
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    Qt = pd.Series(dtype='float64') # initialize a Series to concatenate into
    if pd.notnull(rs.ec):
        year = dt1.year # ONLY do the last year
        this_days = (datetime(year,1,1), datetime(year,12,31))
        print((' ' + str(year) + ' getting ec ' + rn).center(60,'-'))
        rs, qt = rivf.get_ec_data(rs, this_days, temperature=True)
        if rs.got_data:
            qt[qt>100] = np.nan # there are occasional blips
            qt_dict[rn] = qt
        else:
            print(' - no data -')
            
# clean up and organize into a DataFrame
dt00 = datetime(year0,1,1,12)
dt11 = datetime(year1,12,31,12)

t = pd.date_range(start=dt00, end=dt11)
rn_list = list(qt_dict.keys())

new_df = pd.DataFrame(index=t, columns=rn_list)
for rn in rn_list:
    r = qt_dict[rn]
    r = r.sort_index() # the ordering of the historical ec data is off
    r = r.reindex(index=t) # this makes missing dates into nans
    r[r<0] = np.nan # there were a few bad points in willapa 2020
    new_df[rn] = r
    
if not testing:
    # save the result
    new_df.to_pickle(out_dir / ('ALL_temperature.p'))
else:
    print('Testing: new_df not saved.')

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
        new_df[rn].plot(ax=ax)
        ax.set_xlim(dt00, dt11)
        ax.text(.05, .9, rn, transform=ax.transAxes)
        ii += 1

    plt.show()


