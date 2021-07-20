"""
This takes the tidal energy flux through a TEF section and decomposes it into
incident and reflected parts, by tidal frequency.

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
    
# imports
import matplotlib.pyplot as plt
import pickle
import netCDF4 as nc
import pandas as pd
import numpy as np
import utide
from matplotlib.dates import date2num

import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

import plotting_functions as pfun

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun
import flux_fun

# colors
clist = flux_fun.clist

# set input directory
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
dates = '2018.01.01_2018.12.31'
in_dir = in_dir0 / ('extractions_' + dates)

# set input file
# ****************
sect_name = 'ai1'
# ****************
in_fn = in_dir / (sect_name + '.nc')

# get section info
sect_df = tef_fun.get_sect_df(gridname)
x0, x1, y0, y1 = sect_df.loc[sect_name,:]
in_sign = 1 # default, just to get dir_str
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
        
# load extracted fields
ds = nc.Dataset(in_fn)
h = ds['h'][:]
q = ds['q'][:]
DA = ds['DA'][:]
DA0 = ds['DA0'][:]
ot = ds['ocean_time'][:]
zeta = ds['zeta'][:]
NT, NZ, NX = q.shape

# make derived quantities

g = 9.8 # gravity

# mean depth - really we should use the DX or DY-weighted mean
# but it will be very close to this since DX or DY cahnge little
# over a section
mask = ~np.isnan(DA[0,0,:]).data
H = h[mask].mean()

# velocity
omat = np.nan * np.ones(NT)
Q = np.nansum(np.nansum(q, axis=2), axis=1)
A = np.nansum(np.nansum(DA, axis=2), axis=1)
u = Q/A

# normalized velocity [m]
c = np.sqrt(g*H)
un = u / (c/H)

# surface height
eta = zeta[:,mask].mean(axis=1)

# rotate into incident and reflected
th = -45 * np.pi / 180
cth = np.cos(th)
sth = np.sin(th)

R = un*cth + eta*sth

I = eta*cth - un*sth

# harmonic decomposition
def get_harmonics(df, vn, lat):
    t = date2num(df.index.to_pydatetime())
    v = df[vn].to_numpy()
    hm = utide.solve(t, v, v=None,
                 lat=lat,
                 nodal=False,
                 trend=False,
                 method='ols',
                 conf_int='linear',
                 Rayleigh_min=0.95)
    # hm.aux.freq has units cyles/hour
    # so for f = hm.aux.frq[hm.name == 'M2'][0] we get
    # 1/f = 12.420601202671868 (hours per cycle)
    # hm.A is amplitude, hm.g is phase (degrees)
    return hm
    
t_list = []
for t in ot:
    t_list.append(Lfun.modtime_to_datetime(t))
dti = pd.to_datetime(t_list)
dti = dti.tz_localize('UTC')
df = pd.DataFrame(data={'I':I, 'R':R}, index = dti)
df.index.name = 'Date'

lat = (y0 + y1)/2
hm_I = get_harmonics(df, 'I', lat)
hm_R = get_harmonics(df, 'R', lat)

AI = dict()
AR = dict()
gI = dict()
gR = dict()
cons_list = ['M2', 'S2', 'K1', 'O1']
for cons in cons_list:
    if len(cons) == 2: # limit to the simpler constituents
        AI[cons] = hm_I.A[hm_I.name == cons][0]
        AR[cons] = hm_R.A[hm_R.name == cons][0]
        gI[cons] = hm_I.g[hm_I.name == cons][0]
        gR[cons] = hm_R.g[hm_R.name == cons][0]
        print('%s: I=%0.2f (%0.1f deg) R=%0.2f (%0.1f deg) R/I=%0.2f' % (cons, AI[cons], gI[cons], AR[cons], gR[cons], AR[cons]/AI[cons]))

# -----------

# PLOTTING
plt.close('all')
fs = 16 # fontsize
ffs = .7*fs # flux text fontsize
alpha = .5 # transparency for arrows
cf = 'purple' # color for tidal energy flux
cq = 'g' # color for freshwater flux

pfun.start_plot(fs=fs, figsize=(16,8))
fig, axes = plt.subplots(1, 2, squeeze=False)

ax = axes[0,0]
ax.plot(un, eta, '.g', ms=2, alpha=.2)
ax.grid(True)
ax.axhline()
ax.axvline()
ax.set_xlabel('Normalized U')
ax.set_ylabel('Eta')
ax.axis('square')
ax.set_title('%s: %s' % (sect_name, dir_str))

ax = axes[0,1]
ax.plot(R, I, '.g', ms=2, alpha=.2)
ax.grid(True)
ax.axhline()
ax.axvline()
ax.set_xlabel('Reflected')
ax.set_ylabel('Incident')
ax.axis('square')

om = np.linspace(0,2*np.pi,361)
for cons in cons_list:
    ai = AI[cons]
    ar = AR[cons]
    gi = gI[cons] * np.pi / 180
    gr = gR[cons] * np.pi / 180
    if cons[-1] == '1':
        clr = 'b'
    elif cons[-1] == '2':
        clr = 'm'
    ax.plot(ar*np.cos(om+gr), ai*np.cos(om+gi), '-', c=clr, ms=2)
        
plt.show()
pfun.end_plot()