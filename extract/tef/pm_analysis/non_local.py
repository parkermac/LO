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

import Lfun
import zfun
import plotting_functions as pfun

Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtagex', type=str, default='cas6_v3_lo8b')   # e.g. cas6_v3_lo8b
parser.add_argument('-0', '--ds0', type=str, default='2018.01.01')        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='2018.12.31')        # e.g. 2019.07.04
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
sect_list = [item.name.replace('.nc','') for item in in_dir.glob('*.nc')]

for sect_name in ['hc2']:

    # load extracted fields
    ds = nc.Dataset(in_dir / (sect_name + '.nc'))
    q = ds['q'][:]
    s = ds['salt'][:]
    DA = ds['DA'][:]
    ot = ds['ocean_time'][:].data
    td = ot/86400 # time in days
    T = td - td[0] # for plotting

    # precondition fields
    # mask: True for water points on section
    mask = ~np.isnan(q[0,0,:]).data
    q = q[:,:,mask].data
    s = s[:,:,mask].data
    da = DA[:,:,mask].data
    
    # net salt flux
    QS = (q*s).sum(axis=2).sum(axis=1)
    QS_lp = zfun.lowpass(QS, f='godin')
    
    # "tidal" salt flux
    Q = q.sum(axis=2).sum(axis=1)
    Qlp = zfun.lowpass(Q, f='godin')
    Qp = Q - Qlp
    
    sda = s * da
    S = sda.sum(axis=2).sum(axis=1) / da.sum(axis=2).sum(axis=1)
    Slp = zfun.lowpass(S, f='godin')
    Sp = S - Slp
    
    QS_tidal_lp = zfun.lowpass(Qp*Sp, f='godin')
    
    Qlp_Slp = Qlp * Slp
    
    
plt.close('all')
pfun.start_plot(fs=14, figsize=(14,8))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(T, QS_lp/1e3, '-r', T, QS_tidal_lp/1e3, '-b', T, Qlp_Slp/1e3, '-g')
ax.grid(True)
ax.axhline()
ax.set_title(sect_name)

plt.show()
pfun.end_plot()



