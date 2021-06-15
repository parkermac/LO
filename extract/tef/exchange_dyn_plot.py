"""
Plot the exchange flow in a dynamical context.

"""
from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time

import Lfun
import zfun
import plotting_functions as pfun
import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

Ldir = Lfun.Lstart()

# specify section and bulk folder
sect_name = 'ai1'
in_dir = Path('/Users/pm8/Documents/LO_output/extract/cas6_v3t075_lo8/tef/bulk_2018.01.01_2018.12.31')

# get results of bulk_calc
bulk = pickle.load(open(in_dir / (sect_name + '.p'), 'rb'))

# get two-layer time series
gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)

# PLOTTING
fs = 14
pfun.start_plot(fs=fs, figsize=(14,10))

tef_df, in_sign, dir_str = flux_fun.get_two_layer(in_dir, sect_name, gridname)
if tef_df['salt_in'].mean() > tef_df['salt_out'].mean():
    pass
else:
    print('Warning: sign logic breakdown! ' + sect_name)

tef_df['Qe'] = (tef_df['Qin'] - tef_df['Qout'])/2
tef_df['DS'] = (tef_df['salt_in'] - tef_df['salt_out'])
tef_df['Sbar'] = (tef_df['salt_in'] + tef_df['salt_out'])/2

            
fig = plt.figure()

qlim = np.around(1.5*tef_df['Qe'].max(), 0)

# Tranport vs. Time
ax = fig.add_subplot(211)
tef_df[['Qe']].plot(ax=ax, legend=False, color=['r','b'], grid=True)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylim(0, qlim)
ax.set_ylabel('Qe $[m^{3}s^{-1}]$')
ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', color='r',
    transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', color='b',
    transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
# Salinity vs. Time (color by Transport)
ax = fig.add_subplot(212)
tef_df[['DS']].plot(ax=ax, legend=False, color=['r','b'], grid=True)
ax.set_title('Section = ' + sect_name + ': Positive is ' + dir_str)
ax.grid(True)
ax.set_ylabel('Salinity')
ax.text(.03, .95, '(b)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    

    
fig.tight_layout()
plt.show()
pfun.end_plot()