"""
Compare results of LO/extract/tef with those from LiveOcean/x_tef2.

RESULT: they are exactly the same.  Hooray!

"""
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun
import flux_fun
from importlib import reload
reload(flux_fun)

Ldir = Lfun.Lstart()

sect_name = 'ai1'
in_dir_old = Path('/Users/pm8/Documents/LiveOcean_output/tef2/cas6_v3_lo8b_2018.01.01_2018.12.31/bulk')
in_dir = Path('/Users/pm8/Documents/LO_output/extract/cas6_v3_lo8b/tef/bulk_2018.01.01_2018.12.31')

# =================================
# get the DataFrame of all sections
gridname='cas6'
sect_df = tef_fun.get_sect_df(gridname)

# PLOTTING
fs = 14
pfun.start_plot(fs=fs, figsize=(14,10))

tef_df_old, in_sign_old, dir_str, sdir = flux_fun.get_two_layer(in_dir_old, sect_name, gridname, old_style=True)
tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(in_dir, sect_name, gridname)

if in_sign_old != in_sign:
    print('** in_sign error **')
    sys.exit()

# some information about direction
x0, x1, y0, y1 = sect_df.loc[sect_name,:]
if (x0==x1) and (y0!=y1):
    sdir = 'NS'
    if in_sign == 1:
        dir_str = 'Eastward'
    elif in_sign == -1:
        dir_str = 'Westward'
elif (x0!=x1) and (y0==y1):
    sdir = 'EW'
    if in_sign == 1:
        dir_str = 'Northward'
    elif in_sign == -1:
        dir_str = 'Southward'
            
fig = plt.figure()

qlim_p = np.around(1.2*tef_df['Qin'].max(), 0)
qlim_m = np.around(-1.2*tef_df['Qout'].min(), 0)
qlim = np.max([qlim_p, qlim_m])

# Salinity vs. Time (color by Transport)
ax = fig.add_subplot(211)
tef_df[['salt_in','salt_out']].plot(ax=ax, legend=False, color=['r','b'], alpha=.5)
tef_df_old[['salt_in','salt_out']].plot(ax=ax, legend=False, color=['r','b'], linestyle='--')
ax.set_title('Section = ' + sect_name + ': Positive is ' + dir_str)
ax.grid(True)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel('Salinity')
ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
# Tranport vs. Time
ax = fig.add_subplot(212)
tef_df[['Qin','Qout']].plot(ax=ax, legend=False, color=['r','b'], alpha=.5)
tef_df_old[['Qin','Qout']].plot(ax=ax, legend=False, color=['r','b'], linestyle='--')
ax.set_ylim(-qlim, qlim)
ax.grid(True)
ax.set_ylabel('Transport $[m^{3}s^{-1}]$')
ax.text(.03, .95, '(b)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', color='r',
    transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', color='b',
    transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
fig.tight_layout()
plt.show()
pfun.end_plot()