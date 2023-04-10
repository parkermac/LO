"""
Program to plot a map of river tracks used in a specific ROMS grid.

"""

import pandas as pd
import matplotlib.pyplot as plt
import sys

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-ctag', type=str, default='lo_base')
args = parser.parse_args()
ctag = args.ctag

if ctag == 'lo_base':
    gridname = 'cas6'
else:
    print('You need to specify a gridname for this ctag.')
    sys.exit()
    
Ldir = Lfun.Lstart(gridname=gridname)

# get the list of rivers for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')

# find out which rivers have temperature climatologies
riv_dir0 = Ldir['LOo'] / 'pre' / 'river1' / ctag
riv_dir = riv_dir0 / 'Data_historical'
plot_dir = riv_dir0 / 'Data_historical_plots'
Lfun.make_dir(plot_dir)
Ctemp_df = pd.read_pickle(riv_dir / 'CLIM_temp.p')

# PLOTTING

plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(10,13))
fig = plt.figure()

ax = fig.add_subplot(111)
pfun.add_coast(ax, color='gray')
pfun.dar(ax)
ax.set_xlim(-127.5, -121)
ax.set_ylim(42.5, 50.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

for rn in gri_df.index:
    try:
        track_df = pd.read_pickle(riv_dir0 / 'tracks' / (rn + '.p'))
        x = track_df['lon'].to_numpy()
        y = track_df['lat'].to_numpy()
        ax.plot(x, y, '-c', linewidth=2)
        ax.plot(x[-1], y[-1], 'oc', markersize=8)

        if rn in Ctemp_df.columns:
            rcol = 'r'
        else: rcol = 'b'
        
        if rn == 'dosewallips':
            dy = .11
        else:
            dy = .07
        
        ax.text(x[-1]+.06, y[-1]+dy, rn.title(), color=rcol, weight='bold',
            va='center',
            size=.6*fs, rotation=20)
        
    except FileNotFoundError:
        print('no track for ' + rn)
        pass

ax.text(.03,.35,'Rivers for %s' % (gridname), weight='bold', transform=ax.transAxes)
ax.text(.03,.3, 'ctag = %s' % (ctag), weight='bold', transform=ax.transAxes)
ax.text(.03,.25,'RED have Temperature Data', color='r', weight='bold', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(plot_dir / 'river_map.png')
plt.show()
pfun.end_plot()