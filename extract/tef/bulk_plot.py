"""
Plot bulk fluxes as a time series.
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

Ldir = Lfun.Lstart()

in_dir00 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir00)
in_dir0 = in_dir00 / gtagex / 'tef'
ext_name = Lfun.choose_item(in_dir0, tag='bulk', exclude_tag='bulk_plots')
in_dir = in_dir0 / ext_name

sect_list = [item.name for item in in_dir.glob('*.p')]
    
out_dir = in_dir0 / ext_name.replace('bulk', 'bulk_plots')
Lfun.make_dir(out_dir, clean=True)

# =================================
# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df()

testing = True
if testing:
    from importlib import reload
    reload(flux_fun)

if testing:
    sect_list = ['ai1']
else:
    sect_list = list(sect_df.index)

# PLOTTING
fs = 14
pfun.start_plot(fs=fs, figsize=(14,10))

for sect_name in sect_list:

    # ---------------------------------------------------------

    tef_df, in_sign = flux_fun.get_two_layer_all(in_dir, sect_name)
    
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
    tef_df[['salt_in','salt_out']].plot(ax=ax, legend=False, color=['r','b'])
    ax.set_title('Section = ' + sect_name + ': Positive is ' + dir_str)
    ax.grid(True)
    ax.set_xticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('Salinity')
    ax.text(.03, .95, '(a)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        
    # Tranport vs. Time
    ax = fig.add_subplot(212)
    tef_df[['Qin','Qout']].plot(ax=ax, legend=False, color=['r','b'])
    ax.set_ylim(-qlim, qlim)
    ax.grid(True)
    ax.set_ylabel('Transport $[m^{3}s^{-1}]$')
    ax.text(.03, .95, '(b)', va='top', weight='bold', transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))

    if tef_df['salt_in'].mean() > tef_df['salt_out'].mean():
        pass
    else:
        print('Warning: sign logic breakdown! ' + sect_name)
    
    ax.text(.97, .95, 'Inflow', ha='right', va='top', weight='bold', color='r',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
    ax.text(.97, .05, 'Outflow', ha='right', va='bottom', weight='bold', color='b',
        transform=ax.transAxes, size=1.2*fs,
        bbox=dict(facecolor='w', edgecolor='None', alpha=0.5))
        
    fig.tight_layout()
    
    if testing:
        plt.show()
    else:
        plt.savefig(out_dir + sect_name + '.png')
        plt.close()

pfun.end_plot()
