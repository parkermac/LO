"""
Program to plot historical records for rivers.
"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart()

import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime, timedelta
import numpy as np

# Load a dataframe with info for rivers
ri_fn = Ldir['LOo'] / 'pre' / 'river' / 'pnw_all_2021_04' / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

# Load the as-run rivers from forcing files to fill in gaps in ec records
roms_fn = 'cas6_v3_2017.01.01_2020.12.31.p' # copied from LiveOcean_output/river
frc_df = pd.read_pickle(Ldir['LOo'] / 'pre' / 'river' / 'Data_roms' / roms_fn)

# location of historical data to plot
riv_dir = Ldir['LOo'] / 'pre' / 'river' / 'Data_historical'

dt0 = datetime(1980,1,1,12)
dt1 = datetime(2020,12,31,12)
t = pd.date_range(start=dt0, end=dt1)

rn_list = ri_df.index.to_list()

all_df = pd.DataFrame(index=t, columns=rn_list)
for rn in rn_list:
    r = pd.read_pickle(riv_dir / (rn + '.p'))
    r = r.sort_index() # the ordering of the historical ec data is off
    r = r.reindex(index=t) # this makes missing dates into nans
    r[r<0] = np.nan # there were a few bad points in willapa 2020
    all_df[rn] = r
mean_ser = all_df.mean()
mean_ser = mean_ser.sort_values(ascending=False)
all_df = all_df[mean_ser.index] # biggest rivers first

plt.close('all')
ii = 0
for rn in all_df.columns[:9]:
    r = all_df[rn]
    if np.mod(ii,9) == 0:
        fig = plt.figure(figsize=(20,10))
        jj = 1
    ax = fig.add_subplot(3,3,jj)
    r.plot(ax = ax)
    ax.text(.05, .9, '%s: mean = %d $[m^{3}s^{-1}]$' % (rn.title(), int(r.mean())), transform=ax.transAxes, weight='bold')
    ax.set_xlim(dt0,dt1)
    if jj < 7:
        ax.set_xticklabels([])
    jj += 1 # increment panel counter
    ii += 1 # increment figure counter
plt.show()

