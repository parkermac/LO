"""
Compare experiments - in this case for runs that had different tidal forcing.

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

sect_name = 'ai1'
gtagex_list = ['cas6_v3t075_lo8', 'cas6_v3_lo8b', 'cas6_v3t110_lo8']
alpha_dict = dict(zip(gtagex_list, [.3, .7, 1]))
label_dict = dict(zip(gtagex_list, ['75% tide', '100% tide', '110% tide']))

in_dir_dict = dict()
for gtagex in gtagex_list:
    in_dir_dict[gtagex] = Path('/Users/pm8/Documents/LO_output/extract/' + gtagex +
        '/tef/bulk_2018.01.01_2018.12.31')

gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)

# PLOTTING
fs = 14
pfun.start_plot(fs=fs, figsize=(14,10))

tef_df_dict = dict()
in_sign_dict = dict()
dir_str_dict = dict()
for gtagex in gtagex_list:
    tef_df_dict[gtagex], in_sign_dict[gtagex], dir_str_dict[gtagex], sdir = flux_fun.get_two_layer(in_dir_dict[gtagex], sect_name, gridname)
    tef_df_dict[gtagex]['Qin1000'] = tef_df_dict[gtagex]['Qin']/1000
    tef_df_dict[gtagex]['Qout1000'] = tef_df_dict[gtagex]['Qout']/1000
fig = plt.figure()

# Salinity vs. Time (color by Transport)
ax = fig.add_subplot(211)
ii = 0
for gtagex in gtagex_list:
    tef_df_dict[gtagex][['salt_in', 'salt_out']].plot(ax=ax, color=['r', 'b'], alpha=alpha_dict[gtagex], legend=False)
    ax.text(.97, .07*(1 + ii), label_dict[gtagex], ha='right', weight='bold', color='k',
        transform=ax.transAxes, size=1.2*fs, alpha = alpha_dict[gtagex])
    ii += 1
    
ax.set_title('Section = ' + sect_name)
ax.grid(True)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_ylabel('Salinity')
ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
    bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    
# Tranport vs. Time
ax = fig.add_subplot(212)
for gtagex in gtagex_list:
    tef_df_dict[gtagex][['Qin1000','Qout1000']].plot(ax=ax, legend=False, color=['r','b'], alpha=alpha_dict[gtagex])
ax.grid(True)
ax.set_ylabel('Transport $[1000\ m^{3}s^{-1}]$')
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