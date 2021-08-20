"""
Program to plot historical records for rivers.
"""

import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

from lo_tools import Lfun

Ldir = Lfun.Lstart()

gtag = 'cas6_v3'
year0 = 1980
year1 = 2020

# location of historical data to plot
riv_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
all_df = pd.read_pickle(riv_dir / ('ALL_temperature_' + str(year0) + '_' + str(year1) + '.p'))

# Load a dataframe with info for rivers to get
ri_fn = Ldir['LOo'] / 'pre' / 'river' / gtag / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

plt.close('all')

dt0 = datetime(year0,1,1)
dt1 = datetime(year1,12,31)
ii = 0
fig_num = 0
for rn in all_df.columns:
    if pd.notnull(ri_df.loc[rn,'usgs']):
        r = all_df[rn]
        if np.mod(ii,12) == 0:
            fig_num += 1
            fig = plt.figure(figsize=(20,10))
            jj = 1
        ax = fig.add_subplot(3,4,jj)
        r.plot(ax = ax)
        ax.text(.05, .9, '%s $[^{\circ}C]$' % (rn.title()),
            transform=ax.transAxes, weight='bold')
        ax.set_xlim(dt0,dt1)
        ax.set_ylim(0,25)
        if jj <= 8:
            ax.set_xticklabels([])
        jj += 1 # increment panel counter
        ii += 1 # increment figure counter
        if (np.mod(ii,12) == 0) or (rn == all_df.columns[-1]):
            pass
            fig.savefig(riv_dir / ('ALL_temperature_usgs_plot_' + str(fig_num) + '.png'))
            
dt0 = datetime(year1,1,1)
dt1 = datetime(year1,12,31)
ii = 0
fig_num = 0
for rn in all_df.columns:
    if pd.notnull(ri_df.loc[rn,'ec']):
        r = all_df[rn]
        if np.mod(ii,6) == 0:
            fig_num += 1
            fig = plt.figure(figsize=(20,10))
            jj = 1
        ax = fig.add_subplot(2,3,jj)
        r.plot(ax = ax)
        ax.text(.05, .9, '%s $[^{\circ}C]$' % (rn.title()),
            transform=ax.transAxes, weight='bold')
        ax.set_xlim(dt0,dt1)
        ax.set_ylim(0, 25)
        if jj <= 3:
            ax.set_xticklabels([])
        jj += 1 # increment panel counter
        ii += 1 # increment figure counter
        if (np.mod(ii,6) == 0) or (rn == all_df.columns[-1]):
            pass
            fig.savefig(riv_dir / ('ALL_temperature_ec_plot_' + str(fig_num) + '.png'))

plt.show()



