"""
Process TEF extractions, giving transport vs. salinity for:
volume, salt, salinity-squared and other variables.

This is based on LO/extract/tef/process_sections.py but is modified to look
at things like < Net Q * Section averaged S >.


"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# setup
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd

import Lfun
import zfun
import plotting_functions as pfun

Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtagex', type=str, default='cas6_v3_lo8b')          # e.g. cas6_v3_lo8b
parser.add_argument('-0', '--ds0', type=str, default='2018.01.01')        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='2018.12.31')        # e.g. 2019.07.04
parser.add_argument('-sn', '--sect_name', type=str, default='ai1')        # e.g. ai1
args = parser.parse_args()

in_dir00 = Ldir['LOo'] / 'extract'
if len(args.gtagex) == 0:
    gtagex = Lfun.choose_item(in_dir00)
else:
    gtagex = args.gtagex
in_dir0 = in_dir00 / gtagex / 'tef'
if (len(args.ds0)==0) or (len(args.ds1)==0):
    ext_name = Lfun.choose_item(in_dir0, tag='extractions')
else:
    ext_name = 'extractions_' + args.ds0 + '_' + args.ds1
in_dir = in_dir0 / ext_name
#sect_list = [item.name.replace('.nc','') for item in in_dir.glob('*.nc')]

for sect_name in [args.sect_name]:

    # load extracted fields
    ds = nc.Dataset(in_dir / (sect_name + '.nc'))
    q = ds['q'][:]
    s = ds['salt'][:]
    a = ds['DA'][:]
    ot = ds['ocean_time'][:].data
    td = ot/86400 # time in days
    T = td - td[0] # for plotting

    # precondition fields
    # mask: True for water points on section
    mask = ~np.isnan(q[0,0,:]).data
    q = q[:,:,mask].data
    s = s[:,:,mask].data
    a = a[:,:,mask].data
    sa = s*a
    NT, NZ, NX = q.shape
    NC = NZ*NX # number of cells on section
    
    # net salt flux
    F = (q*s).sum(axis=2).sum(axis=1)
    F_lp = zfun.lowpass(F, f='godin')
    
    # standard decomposition
    q_lp = zfun.lowpass(q, f='godin')
    s_lp = zfun.lowpass(s, f='godin')
    sa_lp = zfun.lowpass(sa, f='godin')
    a_lp = zfun.lowpass(a, f='godin')
    Q_lp = q_lp.sum(axis=2).sum(axis=1)
    A_lp = a_lp.sum(axis=2).sum(axis=1)
    S_lp_bar = (sa_lp.sum(axis=2).sum(axis=1)) / A_lp
    
    F0 = Q_lp * S_lp_bar
    
    q_lp_p = q_lp - Q_lp.reshape((NT,1,1))
    s_lp_p = s_lp - S_lp_bar.reshape((NT,1,1))
    
    F1 = ((q_lp_p * s_lp_p).sum(axis=2).sum(axis=1))
    
    F2 = F_lp - F0 - F1
    
    stn_df = pd.DataFrame(index=T)
    stn_df['F_lp'] = F_lp
    stn_df['F0'] = F0
    stn_df['F1'] = F1
    stn_df['F2'] = F2
    stn_df = stn_df/1e3
    
    # "tidal" salt flux
    Q = q.sum(axis=2).sum(axis=1)
    Q_lp = zfun.lowpass(Q, f='godin')
    Q_p = Q - Q_lp
    
    S = sa.sum(axis=2).sum(axis=1) / a.sum(axis=2).sum(axis=1)
    S_lp = zfun.lowpass(S, f='godin')
    S_p = S - S_lp
        
    Qlp_Slp = Q_lp * S_lp
    
    QS_tidal_lp = zfun.lowpass(Q_p*S_p, f='godin')
    
    Q_rest = F_lp - Qlp_Slp - QS_tidal_lp
    
    nolo_df = pd.DataFrame(index=T)
    nolo_df['F_lp'] = F_lp
    nolo_df['Mean'] = Qlp_Slp
    nolo_df['Tidal'] = QS_tidal_lp
    nolo_df['Rest'] = Q_rest
    nolo_df = nolo_df/1e3
    
    
    
plt.close('all')
pfun.start_plot(fs=14, figsize=(14,8))
fig = plt.figure()

ax = fig.add_subplot(211)
stn_df.plot(ax=ax, grid=True)
ax.axhline()
ax.set_title(sect_name)

ax = fig.add_subplot(212)
nolo_df.plot(ax=ax, grid=True)
ax.axhline()

plt.show()
pfun.end_plot()



