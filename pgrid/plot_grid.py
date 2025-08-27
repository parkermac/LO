"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.py.
"""
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pickle
from lo_tools import Lfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',type=str) # e.g. cas6
parser.add_argument('-dmax', default=5, type=int) # max depth for colormap [m]
parser.add_argument('-small', default=False, type=Lfun.boolean_string) # True for laptop size
parser.add_argument('-show_names', default=False, type=Lfun.boolean_string) # True for to show river/traps names
args = parser.parse_args()
zmin = -args.dmax

import gfun
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp

testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfp)

# select grid file
in_fn = gfun.select_file(Gr)

# load the default choices
try:
    dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))
except FileNotFoundError:
    # you could fill this in by hand if you wanted
    dch = {'analytical': False} # hack to make cas6 work

Ldir = Lfun.Lstart(gridname=Gr['gridname'])
# get river info if it exists
do_riv = False
riv_fn0 = Ldir['grid'] / 'river_info.csv' # created by grid_to_LO.py
riv_fn1 = Gr['gdir'] / 'roms_river_info.csv' # created by carve_rivers.py
if riv_fn0.is_file():
    riv_df = pd.read_csv(riv_fn0, index_col='rname')
    do_riv = True
elif riv_fn1.is_file():
    riv_df = pd.read_csv(riv_fn1, index_col='rname')
    do_riv = True


# load the data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = ds.mask_rho.values

lon = ds.lon_rho.values
lat = ds.lat_rho.values

plon, plat = pfun.get_plon_plat(lon,lat)
pad = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-pad, plon[0,-1]+pad, plat[0,0]-pad, plat[-1,0]+pad)

# make a version of z with nans where masked
zm = z.copy()
zm[mask_rho == 0] = np.nan

# PLOTTING
plt.close('all')
if args.small:
    figsize = (8,8)
else:
    figsize = (12,12)
pfun.start_plot(figsize=figsize)

# bathymetry
fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=zmin, vmax=0, cmap='Spectral_r')
# cs = ax.pcolormesh(plon, plat, zm, vmin=-120, vmax=-100, cmap='Spectral_r')
fig.colorbar(cs, ax=ax)
if dch['analytical'] == True:
    pass
else:
    pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(ax_lims)
ax.set_title(in_fn.name)
ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes, bbox=pfun.bbox)
if do_riv:
    gfp.add_riv(riv_df, ds, ax, show_names=args.show_names)
    gfp.add_riv_tracks(Gr, riv_df, ax)

# Also look to see if tiny rivers and wwtp's have been created.
# This is a kludgey step because it depends on having the latest traps names.
traps_list = ['moh20_triv','moh20_wwtp','was24_wwtp'] # from LO/pre/trapsP01 2025.08.27
for traps in traps_list:
    traps_fn = Ldir['grid'] / (traps + '_info.csv')
    if traps_fn.is_file():
        traps_df = pd.read_csv(traps_fn, index_col='rname')
        if 'wwtp' in traps:
            gfp.add_wwtp(traps_df, ds, ax, show_names=args.show_names)
        elif 'triv' in traps:
            gfp.add_triv(traps_df, ds, ax, show_names=args.show_names)

if False:    
    # mask
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tt = ['rho', 'u', 'v']
    sym = dict(zip(['rho', 'u', 'v'],['o','>','^']))
    c = dict(zip(['rho', 'u', 'v'],['b','orange','r']))
    for t in tt:
        x = ds['lon_'+t].values
        y = ds['lat_'+t].values
        m = ds['mask_'+t].values
        ax.plot(x, y, sym[t], c=c[t], alpha=.2, ms=3)
        ax.plot(x[m==1], y[m==1], sym[t], c=c[t], ms=3)
    
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(ax_lims)
    ax.set_title(in_fn.name)
    ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
        
ds.close()

fig.tight_layout()

plt.show()
