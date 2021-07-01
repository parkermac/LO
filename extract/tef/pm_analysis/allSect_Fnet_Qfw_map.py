"""
Plot the mean of tidal energy flux and freshwater transport at all
TEF sections.  It is similar to the supplement figure in MacCready et al. (2021, JGR)
except that it plots freshwater flux instead of volume flux.

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

import Lfun
gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'
#gridname = 'cas6'; tag = 'v3t075'; ex_name = 'lo8'
#gridname = 'cas6'; tag = 'v3t110'; ex_name = 'lo8'
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

import plotting_functions as pfun

pth = Path(__file__).absolute().parent.parent
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun
import flux_fun

# colors
clist = flux_fun.clist

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df(gridname)

# select input directory
in_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef'
dates = '2018.01.01_2018.12.31'
in_fn = in_dir / ('two_layer_mean_' + dates + '.p')
out_fn = in_dir / ('allSect_Fnet_Qfw_map_' + dates + '.png')
# -----------

# PLOTTING
plt.close('all')
fs = 16 # fontsize
ffs = .7*fs # flux text fontsize
alpha = .5 # transparency for arrows
cf = 'purple' # color for tidal energy flux
cq = 'g' # color for freshwater flux

pfun.start_plot(fs=fs, figsize=(10,13))
fig = plt.figure()

# Angles by which to adjust the plotting of transport at each section.
# They are defined -90:90 degrees, positive counter-clockwise, to define the
# direction of an arrow pointing along the local thalweg, relative to the
# normal to the section (positive East or North).
th_dict = {
    'jdf1': -15, 'jdf2': -15, 'jdf3': -5, 'jdf4': 0,
    'sog1': 45, 'sog2': 55, 'sog3': -45, 'sog4': 15, 'sog5': 0,
    'sji1': 0, 'sji2': 20,
    'dp': 0,
    'ai1': -45, 'ai2': 30, 'ai3': 0, 'ai4': -30,
    'wb1': -30, 'wb2': 45, 'wb3': 0, 'wb4': 15,
    'hc1': 10, 'hc2': -45, 'hc3': 10, 'hc4': -50, 'hc5': -45, 'hc6': -30, 'hc7': 0, 'hc8': 25,
    'mb1': 0, 'mb2': -10, 'mb3': 0, 'mb4': 10, 'mb5': -45,
    'tn1': 0, 'tn2': 0, 'tn3': -15,
    'ss1': 45, 'ss2': -30, 'ss3': 0,
}

# get the DataFrame of all sections
sect_df = tef_fun.get_sect_df(gridname)
sect_list = list(sect_df.index)

df = pd.read_pickle(in_fn)
df['Fnet_MW'] = df['Fnet']/1e6 # tidal energy flux in MegaWatts
#
#
# # axis limits
x0 = -125.5; x1 = -122; y0 = 47; y1 = 50.5 # Salish Sea
x00 = -123.3; x11 = -122.2; y00 = 47; y11 = 48.5 # Puget Sound
aaS = [x0, x1, y0, y1]
aaP = [x00, x11, y00, y11]

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
#
for sn in df.index:
    sx = df.loc[sn,'lon']
    sy = df.loc[sn,'lat']
    sdir = df.loc[sn,'sdir']
    landward = df.loc[sn,'in_sign']
    F = df.loc[sn,'Fnet_MW']
    Q = df.loc[sn,'Qfw'] # could also easily substitute Qnet to look at volume flux
    logF = np.log10(np.abs(F) + 1)
    logQ = np.log10(np.abs(Q) + 1)
    sgnF = np.sign(landward*F)
    sgnQ = np.sign(landward*Q)
    scl = 40 # a vector of length 1 will be 1/scl of the length of the y-axis
    if sdir == 'EW':
        om = 90 + th_dict[sn]
    elif sdir == 'NS':
        om = th_dict[sn]
    sth = np.sin(np.pi*om/180)
    cth = np.cos(np.pi*om/180)
    vf0 = cth*sgnF*logF; vf1 = sth*sgnF*logF
    vq0 = cth*sgnQ*logQ; vq1 = sth*sgnQ*logQ

    ax1.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
        linewidths=1, color=cf, alpha=alpha)
    if sx < x00 or sy > y11:
        ax1.text(sx, sy+.04, str(int(np.abs(F))), ha='center', va='center', size=ffs,
            weight='bold', color=cf, alpha=1)

    if sx > x00 and sy < y11:
        ax2.quiver(sx,sy, vf0, vf1, scale=scl, scale_units='height',
            linewidths=1, color=cf, alpha=alpha)
        ax2.text(sx, sy+.02, str(int(np.abs(F))), ha='center', va='center', size=ffs,
        weight='bold', color=cf, alpha=1)

    ax3.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
        linewidths=1, color=cq, alpha=alpha)
    if sx < x00 or sy > y11:
        ax3.text(sx, sy+.04, str(int(np.abs(Q))), ha='center', va='center', size=ffs,
            weight='bold', color=cq, alpha=1)

    if sx > x00 and sy < y11:
        ax4.quiver(sx,sy, vq0, vq1, scale=scl, scale_units='height',
            linewidths=1, color=cq, alpha=alpha)
        ax4.text(sx, sy+.02, str(int(np.abs(Q))), ha='center', va='center', size=ffs,
        weight='bold', color=cq, alpha=1)

pfun.add_coast(ax1, color='gray')
ax1.axis(aaS)
pfun.dar(ax1)
ax1.set_xticklabels([])
ax1.tick_params(labelsize=.8*fs)
ax1.text(.95,.9,'(a)', size=fs, transform=ax1.transAxes, ha='right', weight='bold')
ax1.set_ylabel('Latitude', size=.8*fs)
ax1.set_xticks([-125, -124, -123, -122])
ax1.text(.05,.05,'Tidal energy Flux $[MW]$', size=.8*fs, transform=ax1.transAxes,
    weight='bold', color=cf)

pfun.add_coast(ax2, color='gray')
ax2.axis(aaP)
pfun.dar(ax2)
ax2.set_xticklabels([])
ax2.set_yticks([47, 48])
ax2.tick_params(labelsize=.8*fs)
ax2.text(.95,.9,'(b)', size=fs, transform=ax2.transAxes, ha='right', weight='bold')
ax2.set_xticks([-123, -122.5])

pfun.add_coast(ax3, color='gray')
ax3.axis(aaS)
pfun.dar(ax3)
ax3.tick_params(labelsize=.8*fs)
ax3.text(.95,.9,'(c)', size=fs, transform=ax3.transAxes, ha='right', weight='bold')
ax3.set_xlabel('Longitude', size=.8*fs)
ax3.set_ylabel('Latitude', size=.8*fs)
ax3.set_xticks([-125, -124, -123, -122])
ax3.set_xticklabels([-125, -124, -123, -122])
ax3.text(.05,.05,'Freshwater Flux $[m^{3}s^{-1}]$', size=.8*fs, transform=ax3.transAxes,
    weight='bold', color=cq)

pfun.add_coast(ax4, color='gray')
ax4.axis(aaP)
pfun.dar(ax4)
ax4.set_yticks([47, 48])
ax4.tick_params(labelsize=.8*fs)
ax4.text(.95,.9,'(d)', size=fs, transform=ax4.transAxes, ha='right', weight='bold')
ax4.set_xlabel('Longitude', size=.8*fs)
ax4.set_xticks([-123, -122.5])
ax4.set_xticklabels([-123, -122.5])

fig.suptitle(Ldir['gtagex'] + ': ' + in_fn.name.replace('.p','').replace('two_layer_mean_',''))
fig.tight_layout()

fig.savefig(out_fn)
plt.show()
pfun.end_plot()