"""
Program to gather historical temperature records for rivers.
"""

from lo_tools import Lfun, zfun
from lo_tools import river_functions as rivf
from importlib import reload
reload(rivf)

from datetime import datetime, timedelta
import pandas as pd
import numpy as np

Ldir = Lfun.Lstart()

# set the gtag for this case (could make it a command line argument)
gtag = 'cas6_v3'

# Load a dataframe with info for rivers to get
ri_fn = Ldir['LOo'] / 'pre' / 'river' / gtag / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

# location for output
out_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
Lfun.make_dir(out_dir)

# set time range
year0 = 1980
year1 = 2020

testing = True

if testing:
    ri_df = ri_df.loc[['skagit','fraser']]
    year0 = 2019
    year1 = 2020
        
dt0 = datetime(year0,1,1)
dt1 = datetime(year1,12,31)
days = (dt0, dt1)

qt_dict = dict()
# get USGS river data
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.usgs):
        print('\n'+(' getting usgs ' + rn).center(60,'-'))
        rs, qt = rivf.get_usgs_data(rs, days, temperature=True)
        if rs['got_data']:
            qt_dict[rn] = qt
        else:
            print(' - no data -')

# get ec data, just for most recent year
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    Qt = pd.Series(dtype='float64') # initialize a Series to concatenate into
    if pd.notnull(rs.ec):
        year = dt1.year
        this_days = (datetime(year,1,1), datetime(year,12,31))
        print((' ' + str(year) + ' getting ec ' + rn).center(60,'-'))
        rs, qt = rivf.get_ec_data(rs, this_days, temperature=True)
        if rs['got_data']:
            qt[qt>100] = np.nan # there are occasional blips
            qt_dict[rn] = qt
        else:
            print(' - no data -')
                
            
# clean up and organize into a DataFrame
dt00 = datetime(year0,1,1,12)
dt11 = datetime(year1,12,31,12)

t = pd.date_range(start=dt00, end=dt11)
rn_list = list(qt_dict.keys()) #ri_df.index.to_list()

all_df = pd.DataFrame(index=t, columns=rn_list)
for rn in rn_list:
    r = qt_dict[rn]
    r = r.reindex(index=t) # this makes missing dates into nans
    r[r<0] = np.nan # there were a few bad points in willapa 2020
    all_df[rn] = r
all_df.to_pickle(out_dir / ('ALL_temperature_' + str(year0) + '_' + str(year1) + '.p'))

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


