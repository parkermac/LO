"""
Program to plot historical temperature records for rivers.
"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import zfun
Ldir = Lfun.Lstart()

import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
import numpy as np

# Load a dataframe with info for rivers
ri_fn = Ldir['LOo'] / 'pre' / 'river' / 'pnw_all_2021_04' / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

# location of historical data to plot
riv_dir = Ldir['LOo'] / 'pre' / 'river' / 'Data_historical'

# set separate time limits (mirrors what we used in make_historical_temperature.py)
usgs_dt0 = datetime(1980,1,1,12)
usgs_dt1 = datetime(2020,12,31,12)
ec_dt0 = datetime(2020,1,1,12)
ec_dt1 = datetime(2020,12,31,12)

# usgs stations
t = pd.date_range(start=usgs_dt0, end=usgs_dt1)
all_usgs_df = pd.DataFrame(index=t)
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.usgs):
        r_fn = riv_dir / (rn + '_temperature.p')
        if r_fn.is_file():
            r = pd.read_pickle(riv_dir / (rn + '_temperature.p'))
            r = r.reindex(index=t) # this makes missing dates into nans
            all_usgs_df[rn] = r
all_usgs_df.to_pickle(riv_dir / 'ALL_temperature_data_usgs_df.p')
            
# ec stations
dt0 = datetime(2020,1,1,12)
dt1 = datetime(2020,12,31,12)
t = pd.date_range(start=ec_dt0, end=ec_dt1)
all_ec_df = pd.DataFrame(index=t)
for rn in ri_df.index:
    rs = ri_df.loc[rn].copy() # a series with info for this river
    if pd.notnull(rs.ec):
        r_fn = riv_dir / (rn + '_temperature.p')
        if r_fn.is_file():
            r = pd.read_pickle(riv_dir / (rn + '_temperature.p'))
            r = r.reindex(index=t) # this makes missing dates into nans
            all_ec_df[rn] = r
all_ec_df.to_pickle(riv_dir / 'ALL_temperature_data_ec_df.p')

# PLOTTING
plt.close('all')

# usgs rivers
NP = len(all_usgs_df.columns)
NR, NC = zfun.get_rc(NP)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
ii = 0
for rn in all_usgs_df.columns:
    ir, ic = zfun.get_irc(ii, NC)
    ax = axes[ir, ic]
    r = all_usgs_df[rn]
    r.plot(ax = ax)
    ax.text(.05, .9, '%s Temperature [degC]' % (rn.title()), transform=ax.transAxes, weight='bold')
    ax.set_xlim(usgs_dt0,usgs_dt1)
    ax.set_ylim(0,25)
    if ir < NR-1:
        ax.set_xticklabels([])
    ii += 1 # increment panel counter
fig.tight_layout()
fig.savefig(riv_dir / 'ALL_temperature_plot_usgs.png')

# ec rivers
NP = len(all_ec_df.columns)
NR, NC = zfun.get_rc(NP)
fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)
ii = 0
for rn in all_ec_df.columns:
    ir, ic = zfun.get_irc(ii, NC)
    ax = axes[ir, ic]
    r = all_ec_df[rn]
    r.plot(ax = ax)
    ax.text(.05, .9, '%s Temperature [degC]' % (rn.title()), transform=ax.transAxes, weight='bold')
    ax.set_xlim(ec_dt0,ec_dt1)
    ax.set_ylim(0,25)
    if ir < NR-1:
        ax.set_xticklabels([])
    ii += 1 # increment panel counter
fig.tight_layout()
fig.savefig(riv_dir / 'ALL_temperature_plot_ec.png')

plt.show()

