"""
Program to plot a map of river tracks used in a specific ROMS grid.

"""
import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart(gridname='cas6', tag='v3')
import plotting_functions as pfun

import pandas as pd
import matplotlib.pyplot as plt

# get the list of rivers for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')

# find out which rivers have temperature climatologies
ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
year0 = 1980
year1 = 2020
Ctemp_df = pd.read_pickle(ri_dir / 'Data_historical' /
    ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p'))

# PLOTTING

plt.close('all')
fs=14
plt.rc('font', size=fs)

fig = plt.figure(figsize=(10,13))
ax = fig.add_subplot(111)
pfun.add_coast(ax, color='gray')
pfun.dar(ax)
ax.set_xlim(-127.5, -121)
ax.set_ylim(42.5, 50.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

for rn in gri_df.index:
    try:
        track_df = pd.read_pickle(ri_dir / 'tracks' / (rn + '.p'))
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

ax.text(.03,.35,'LiveOcean Rivers', weight='bold', transform=ax.transAxes)
ax.text(.03,.3, Ldir['gtag'], weight='bold', transform=ax.transAxes)
ax.text(.03,.25,'RED have Temperature Data', color='r', weight='bold', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(ri_dir / 'river_map.png')
plt.show()
plt.rcdefaults()