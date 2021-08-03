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
import Lfun
import zfun
import plotting_functions as pfun

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import flux_fun
  
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtagex', type=str, default='cas6_v3_lo8b')          # e.g. cas6_v3_lo8b
parser.add_argument('-0', '--ds0', type=str, default='2018.01.01')        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='2018.12.31')        # e.g. 2019.07.04
parser.add_argument('-sn', '--sect_name', type=str, default='ai1')        # e.g. ai1
args = parser.parse_args()

Ldir = Lfun.Lstart()
sect_name = args.sect_name
date_range = args.ds0 + '_' + args.ds1
in_dir0 = Ldir['LOo'] / 'extract' / args.gtagex / 'tef'

# get TEF time_series
tef_in_dir = in_dir0 / ('bulk_' + date_range)
# get two-layer time series
tef_df, in_sign, dir_str, sdir = flux_fun.get_two_layer(tef_in_dir, sect_name, 'cas6')
# make derived variables
tef_df['Qe'] = ((tef_df['Qin'] - tef_df['Qout'])/2) / 1e3
tef_df['DS'] = tef_df['salt_in'] - tef_df['salt_out']
tef_df['QeDS'] = tef_df['Qe'] * tef_df['DS']
tef_df['Sbar'] = (tef_df['salt_in'] + tef_df['salt_out'])/2
tef_df['QrSbar'] = tef_df['qnet'] * tef_df['Sbar'] / 1e3
tef_df['Qprism'] = (tef_df['qabs']/2)
# use Freshwater Flux as an alternate way to calculate Qr
Socn = 34
tef_df['Qfw'] = (tef_df['Qin']*(Socn-tef_df['salt_in']) + tef_df['Qout']*(Socn-tef_df['salt_out']))/Socn
# also calculate the related Saltwater Flux
tef_df['Qsw'] = (tef_df['Qin']*tef_df['salt_in'] + tef_df['Qout']*tef_df['salt_out'])/Socn

# Load and process original hourly extracted fields

ext_in_dir = in_dir0 / ('extractions_' + date_range)
ds = nc.Dataset(ext_in_dir / (sect_name + '.nc'))
q = ds['q'][:] # packed as (time, z, x)
s = ds['salt'][:]
da = ds['DA'][:]
ot = ds['ocean_time'][:].data
td = ot/86400 # time in days
T = td - td[0] # for plotting
# precondition fields
# mask: True for water points on section (x)
mask = ~np.isnan(q[0,0,:]).data
q = in_sign * q[:,:,mask].data/1e3
s = s[:,:,mask].data
da = da[:,:,mask].data
u = q/da
sda = s*da
NT, NZ, NX = q.shape
NC = NZ*NX # number of  z,x cells on section

# net salt flux
Fraw = (q*s).sum(axis=2).sum(axis=1)
# - subsample and trim in the same way as bulk_calc.py to match tef_df
pad = 36
F = zfun.lowpass(Fraw, f='godin')[pad:-pad+1:24]

# Standard decomposition

# Use the notation of MacCready and Banas (2011) Treatise chapter.
q0 = zfun.lowpass(q, f='godin')[pad:-pad+1:24]
s0 = zfun.lowpass(s, f='godin')[pad:-pad+1:24]
sda0 = zfun.lowpass(sda, f='godin')[pad:-pad+1:24]
da0 = zfun.lowpass(da, f='godin')[pad:-pad+1:24]
Q0 = q0.sum(axis=2).sum(axis=1)
A0 = da0.sum(axis=2).sum(axis=1)
S0 = (sda0.sum(axis=2).sum(axis=1)) / A0
U0 = Q0/A0

F0 = Q0 * S0

NT0 = len(F0)
u1 = q0/da0 - U0.reshape((NT0,1,1))
s1 = sda0/da0 - S0.reshape((NT0,1,1))

F1 = ((u1 * s1 * da0).sum(axis=2).sum(axis=1))

F2 = F - F0 - F1

stnd_df = pd.DataFrame(index=tef_df.index)
stnd_df['F'] = F
stnd_df['F0'] = F0
stnd_df['F1'] = F1
stnd_df['F2'] = F2
# add a column from the TEF version, for plotting
stnd_df['QeDS'] = tef_df['QeDS']

# # "non-local" salt flux
Q = q.sum(axis=2).sum(axis=1)
A = da.sum(axis=2).sum(axis=1)
S = sda.sum(axis=2).sum(axis=1)/A

Qlp = zfun.lowpass(Q, f='godin')
Slp = zfun.lowpass(S, f='godin')

Qp = Q - Qlp
Sp = S - Slp

F00 = (Qlp * Slp)[pad:-pad+1:24]
F11 = zfun.lowpass((Qp * Sp), f='godin')[pad:-pad+1:24]

stnd_df['F00'] = F00
stnd_df['Fnonlocal'] = F11

# Q = q.sum(axis=2).sum(axis=1)
# Q_lp = zfun.lowpass(Q, f='godin')
# Q_p = Q - Q_lp
#
# S = sa.sum(axis=2).sum(axis=1) / a.sum(axis=2).sum(axis=1)
# S_lp = zfun.lowpass(S, f='godin')
# S_p = S - S_lp
#
# Qlp_Slp = Q_lp * S_lp
#
# QS_tidal_lp = zfun.lowpass(Q_p*S_p, f='godin')
#
# Q_rest = F_lp - Qlp_Slp - QS_tidal_lp
#
# nolo_df = pd.DataFrame(index=tef_df.index)
# nolo_df['F_lp'] = F_lp
# nolo_df['Mean'] = Qlp_Slp
# nolo_df['Tidal'] = QS_tidal_lp
# nolo_df['Rest'] = Q_rest
# nolo_df = nolo_df/1e3

# PLOTTING
plt.close('all')
pfun.start_plot(fs=14, figsize=(14,8))
fig = plt.figure()

ax = fig.add_subplot(211)
stnd_df[['F1','F2','QeDS']].plot(ax=ax, grid=True)
ax.axhline(c='k')
ax.set_title(sect_name)
ax.set_xticklabels([])

ax = fig.add_subplot(212)
stnd_df[['F1','F2', 'Fnonlocal']].plot(ax=ax, grid=True)
ax.axhline(c='k')

plt.show()
pfun.end_plot()



