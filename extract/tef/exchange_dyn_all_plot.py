"""
Plot the exchange flow in a dynamical context, using the time mean
of all sections

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
import pandas as pd
from time import time

import Lfun
import zfun
import plotting_functions as pfun
import tef_fun
import flux_fun

Ldir = Lfun.Lstart()

gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)

ii = 0
for gtagex in ['cas6_v3_lo8b', 'cas6_v3t075_lo8']:
    in_dir = Path('/Users/pm8/Documents/LO_output/extract/'+gtagex+'/tef/bulk_2018.01.01_2018.12.31')
    tt0 = time()
    df = pd.DataFrame(index=sect_list)
    for sect_name in sect_list:
        gridname = 'cas6'
        tef_df, in_sign, dir_str = flux_fun.get_two_layer(in_dir, sect_name, gridname)
        # make derived variables
        df.loc[sect_name, 'Qe'] = ((tef_df['Qin'] - tef_df['Qout']).mean()/2)/1000
        df.loc[sect_name, 'DS'] = (tef_df['salt_in'] - tef_df['salt_out']).mean()
        df.loc[sect_name, 'Sbar'] = (tef_df['salt_in'] + tef_df['salt_out']).mean()/2
        df.loc[sect_name, 'Qprism'] =( tef_df['qabs'].mean()/2)/1000
        df.loc[sect_name, 'SSH'] = tef_df['ssh'].mean()
        df.loc[sect_name, 'Fnet'] = tef_df['fnet'].mean()
        df.loc[sect_name, 'Qnet'] = tef_df['qnet'].mean()
    if ii == 0:
        df0 = df.copy()
    elif ii == 1:
        df1 = df.copy()
    print('Total time to fill DataFrame = %0.1f sec' % (time()-tt0))
    ii += 1

# PLOTTING
c0 = 'darkred'
c1 = 'dodgerblue'
alpha = 0.5
loglog = True

plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(12,12))

ax = df0.plot(x='Qprism', y='Qe', linestyle='None', marker='o',
    color=c0, label='Original', alpha=alpha, loglog=loglog)
df1.plot(x='Qprism', y='Qe', linestyle='None', marker='o',
    color=c1, ax=ax, label='75% tide', alpha=alpha, loglog=loglog)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$Q_{E}\ [10^{3}\ m^{3}s^{-1}]$')

for sect_name in sect_list:
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'Qe'], sect_name, fontsize=.7*(fs), color=c0)
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Qe'], df1.loc[sect_name,'Qe']], '-', c='gray')

plt.show()
pfun.end_plot()

