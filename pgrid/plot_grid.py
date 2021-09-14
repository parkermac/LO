"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.py.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',
        type=str)
args = parser.parse_args()

import gfun
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp

import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfp)

# select grid file
in_fn = gfun.select_file(Gr)

# get river info if it exists
do_riv = False
ri_fn = Gr['gdir'] / 'roms_river_info.csv'
if ri_fn.is_file():
    rri_df = pd.read_csv(ri_fn, index_col='rname')
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
pfun.start_plot(figsize=(12,12))
fig = plt.figure()

ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=-150, vmax=10, cmap='Spectral_r')
fig.colorbar(cs, ax=ax)
pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(ax_lims)
ax.set_title(in_fn.name)
ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes)

if do_riv:
    gfp.add_river_tracks(Gr, ds, ax)
        
ds.close()

plt.show()
