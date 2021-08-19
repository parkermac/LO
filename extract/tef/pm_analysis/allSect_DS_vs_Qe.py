"""
Plot DS vs. Qe time series as a little scribble for each section.

"""
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from time import time

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun
import flux_fun

testing = True

gtagex = 'cas6_v3_lo8b'
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)
if testing:
    sect_list = ['ai1', 'ai4', 'mb3', 'tn2']

# specify bulk folder
in_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / 'bulk_2018.01.01_2018.12.31'

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(24,8))

fig = plt.figure()

ax1 = fig.add_subplot(131)
ax1.grid(True)
ax1.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax1.set_ylabel(r'$\Delta S$')

ax2 = fig.add_subplot(132)
ax2.grid(True)
ax2.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax2.set_ylabel(r'$Q_{E}\ [10^{3}\ m^{3}s^{-1}]$')

ax3 = fig.add_subplot(133)
ax3.grid(True)
ax3.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax3.set_ylabel(r'$Q_{E} \Delta S\ [10^{3}\ m^{3}s^{-1}]$')

for sect_name in sect_list:
    # get two-layer time series
    tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(in_dir, sect_name, 'cas6')
    # make derived variables
    tef_df['Qe'] = ((tef_df['Qin'] - tef_df['Qout'])/2)/1000
    tef_df['DS'] = tef_df['salt_in'] - tef_df['salt_out']
    tef_df['QeDS'] = tef_df['Qe'] * tef_df['DS']
    tef_df['Sbar'] = (tef_df['salt_in'] + tef_df['salt_out'])/2
    tef_df['Qprism'] = (tef_df['qabs']/2)/1000
    # use Freshwater Flux as an alternate way to calculate Qr
    Socn = 34
    tef_df['Qfw'] = (tef_df['Qin']*(Socn-tef_df['salt_in']) + tef_df['Qout']*(Socn-tef_df['salt_out']))/Socn
    # also calculate the related Saltwater Flux
    tef_df['Qsw'] = (tef_df['Qin']*tef_df['salt_in'] + tef_df['Qout']*tef_df['salt_out'])/Socn
    # drop times with negative DS
    tef_df.loc[tef_df['DS'] <= 0, :] = np.nan
    # and plot this section
    ax1.plot(tef_df['Qprism'].to_numpy(),tef_df['DS'].to_numpy(),'-', label=sect_name, alpha=.8)
    ax2.plot(tef_df['Qprism'].to_numpy(),tef_df['Qe'].to_numpy(),'-', label=sect_name, alpha=.8)
    ax3.plot(tef_df['Qprism'].to_numpy(),(tef_df['Qe']*tef_df['DS']).to_numpy(),'-', label=sect_name, alpha=.8)

    
ax1.legend()
fig.tight_layout()
plt.show()
pfun.end_plot()
