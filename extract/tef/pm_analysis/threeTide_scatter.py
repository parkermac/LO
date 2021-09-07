"""
Plot the exchange flow in a dynamical context, using the time mean
of all sections

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

Ldir = Lfun.Lstart()

gridname = 'cas6'
sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)

ii = 0
for gtagex in ['cas6_v3_lo8b', 'cas6_v3t075_lo8', 'cas6_v3t110_lo8']:
    in_fn = Path('/Users/pm8/Documents/LO_output/extract/' + gtagex +
        '/tef/two_layer_mean_2018.08.01_2018.12.31.p')
    df = pd.read_pickle(in_fn)
    if ii == 0:
        df0 = df.copy()
    elif ii == 1:
        df1 = df.copy()
    elif ii == 2:
        df2 = df.copy()
    ii += 1

# PLOTTING
c0 = 'g'
c1 = 'dodgerblue'
c2 = 'darkred'
alpha = 0.4
loglog = True
logx = True
ms = 8
fs = 14

plt.close('all')
pfun.start_plot(fs=fs, figsize=(10,10))

# Qe vs. Qprism
ax = df0.plot(x='Qprism', y='Qe', linestyle='None', marker='o',
    color=c0, label='Original', alpha=alpha, loglog=loglog, markersize=ms)
df1.plot(x='Qprism', y='Qe', linestyle='None', marker='o',
    color=c1, ax=ax, label='75% tide', alpha=alpha, loglog=loglog, markersize=ms)
df2.plot(x='Qprism', y='Qe', linestyle='None', marker='o',
    color=c2, ax=ax, label='110% tide', alpha=alpha, loglog=loglog, markersize=ms)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$Q_{E}\ [10^{3}\ m^{3}s^{-1}]$')
for sect_name in sect_list:
    # add section names
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'Qe'], sect_name, fontsize=.7*(fs),
        color='k', ha='center', va='center')
    # add lines connecting the experiments for each section
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Qe'], df1.loc[sect_name,'Qe']], '-', c='gray')
    ax.plot([df0.loc[sect_name,'Qprism'], df2.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Qe'], df2.loc[sect_name,'Qe']], '-', c='gray')
            
# Qfw vs. Qprism
ax = df0.plot(x='Qprism', y='Qfw', linestyle='None', marker='o',
    color=c0, label='Original', alpha=alpha, logx=logx, markersize=ms)
df1.plot(x='Qprism', y='Qfw', linestyle='None', marker='o',
    color=c1, ax=ax, label='75% tide', alpha=alpha, logx=logx, markersize=ms)
df2.plot(x='Qprism', y='Qfw', linestyle='None', marker='o',
    color=c2, ax=ax, label='110% tide', alpha=alpha, logx=logx, markersize=ms)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$Q_{FW}\ [10^{3}\ m^{3}s^{-1}]$')
for sect_name in sect_list:
    # add section names
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'Qfw'], sect_name, fontsize=.7*(fs),
        color='k', ha='center', va='center')
    # add lines connecting the experiments for each section
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Qfw'], df1.loc[sect_name,'Qfw']], '-', c='gray')
    ax.plot([df0.loc[sect_name,'Qprism'], df2.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Qfw'], df2.loc[sect_name,'Qfw']], '-', c='gray')
            
# Sbar vs. Qprism
ax = df0.plot(x='Qprism', y='Sbar', linestyle='None', marker='o', logx=logx,
    color=c0, label='Original', alpha=alpha, markersize=ms)
df1.plot(x='Qprism', y='Sbar', linestyle='None', marker='o', logx=logx,
    color=c1, ax=ax, label='75% tide', alpha=alpha, markersize=ms)
df2.plot(x='Qprism', y='Sbar', linestyle='None', marker='o', logx=logx,
    color=c2, ax=ax, label='110% tide', alpha=alpha, markersize=ms)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$Sbar$')
for sect_name in sect_list:
    # add section names
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'Sbar'], sect_name, fontsize=.7*(fs),
        color='k', ha='center', va='center')
    # add lines connecting the experiments for each section
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Sbar'], df1.loc[sect_name,'Sbar']], '-', c='gray')
    ax.plot([df0.loc[sect_name,'Qprism'], df2.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'Sbar'], df2.loc[sect_name,'Sbar']], '-', c='gray')
            
# DS vs. Qprism
ax = df0.plot(x='Qprism', y='DS', linestyle='None', marker='o', logx=logx,
    color=c0, label='Original', alpha=alpha, markersize=ms)
df1.plot(x='Qprism', y='DS', linestyle='None', marker='o', logx=logx,
    color=c1, ax=ax, label='75% tide', alpha=alpha, markersize=ms)
df2.plot(x='Qprism', y='DS', linestyle='None', marker='o', logx=logx,
    color=c2, ax=ax, label='110% tide', alpha=alpha, markersize=ms)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$\Delta S$')
for sect_name in sect_list:
    # add section names
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'DS'], sect_name, fontsize=.7*(fs),
        color='k', ha='center', va='center')
    # add lines connecting the experiments for each section
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'DS'], df1.loc[sect_name,'DS']], '-', c='gray')
    ax.plot([df0.loc[sect_name,'Qprism'], df2.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'DS'], df2.loc[sect_name,'DS']], '-', c='gray')
            
# Qin*DS vs. Qprism
df0['QinDS'] = df0['Qin'] * df0['DS']
df1['QinDS'] = df1['Qin'] * df1['DS']
df2['QinDS'] = df2['Qin'] * df2['DS']
ax = df0.plot(x='Qprism', y='QinDS', linestyle='None', marker='o', loglog=loglog,
    color=c0, label='Original', alpha=alpha, markersize=ms)
df1.plot(x='Qprism', y='QinDS', linestyle='None', marker='o', loglog=loglog,
    color=c1, ax=ax, label='75% tide', alpha=alpha, markersize=ms)
df2.plot(x='Qprism', y='QinDS', linestyle='None', marker='o', loglog=loglog,
    color=c2, ax=ax, label='110% tide', alpha=alpha, markersize=ms)
ax.set_xlabel(r'$Q_{prism}\ [10^{3}\ m^{3}s^{-1}]$')
ax.set_ylabel(r'$Qin\Delta S$')
for sect_name in sect_list:
    # add section names
    ax.text(df0.loc[sect_name,'Qprism'], df0.loc[sect_name,'QinDS'], sect_name, fontsize=.7*(fs),
        color='k', ha='center', va='center')
    # add lines connecting the experiments for each section
    ax.plot([df0.loc[sect_name,'Qprism'], df1.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'QinDS'], df1.loc[sect_name,'QinDS']], '-', c='gray')
    ax.plot([df0.loc[sect_name,'Qprism'], df2.loc[sect_name,'Qprism']],
            [df0.loc[sect_name,'QinDS'], df2.loc[sect_name,'QinDS']], '-', c='gray')


plt.show()
pfun.end_plot()

