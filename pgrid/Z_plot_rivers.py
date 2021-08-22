# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:10:34 2016

@author: PM5

Code to plot river tracks and names, and then associate the names with
rivers that have temperature climatology data.

"""

import os
import sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))

from importlib import reload
import Lfun
Ldir = Lfun.Lstart()

sys.path.append(os.path.abspath('../../LiveOcean/plotting'))
import pfun

sys.path.append(os.path.abspath('../../LiveOcean/forcing/riv2/'))
import river_functions as rivfun

import gfun
G = gfun.gstart()

import pandas as pd
import matplotlib.pyplot as plt

# get river info
ri_fn = G['ri_dir'] + 'river_info.csv'
df = pd.read_csv(ri_fn, index_col='rname')
dir0 = Ldir['data'] + 'rivers/'
clim_dir = dir0 + 'Data_T_clim/'
c_list_raw = os.listdir(clim_dir)
c_list = []
for m in c_list_raw:
    if '.csv' in m:
        c_list.append(m)

# associate rivers with ones that have temperature climatology data
df = rivfun.get_tc_rn(df)

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

for rn in df.index:
    try:
        fn_tr = G['ri_dir'] + 'tracks/' + rn + '.csv'
        df_tr = pd.read_csv(fn_tr, index_col='ind')
        x = df_tr['lon'].to_numpy()
        y = df_tr['lat'].to_numpy()
        ax.plot(x, y, '-c', linewidth=2)
        ax.plot(x[-1], y[-1], 'oc', markersize=8)

        if (rn + '.csv') in c_list:
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
        pass

ax.text(.03,.3,'LiveOcean Rivers', weight='bold', transform=ax.transAxes)
ax.text(.03,.25,'RED have Temp. Data', color='r', weight='bold', transform=ax.transAxes)
fig.tight_layout()
plt.show()
plt.rcdefaults()