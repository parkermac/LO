"""
Code to make an exploratory plot of the output of a box extraction.
"""
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys
from cmocean import cm
import pandas as pd

Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'box'
box_name = Lfun.choose_item(in_dir, tag='.nc', exclude_tag='', itext='** Choose box extraction from list **')
box_fn = in_dir / box_name

# gather fields
ds = xr.open_dataset(box_fn)

# get time
ot = ds.ocean_time.values
ot_dt = pd.to_datetime(ot)
NT = len(ot_dt)
# choose a time to get
print('Time range = ')
print('0 = %s UTC' % (ot_dt[0].strftime('%Y.%m.%d %H:%M:%S')))
print('%d = %s UTC' % (NT-1, ot_dt[-1].strftime('%Y.%m.%d %H:%M:%S')))
my_choice = input('-- Input time index to plot -- (return=0) ')
if len(my_choice)==0:
    my_choice = 0
if int(my_choice) not in range(NT):
    print('Error: time index out of range.')
    sys.exit()
else:
    nt = int(my_choice)

lon = ds.lon_rho.values
lat = ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon, lat)

plot_uv = False
if 'u' in ds.data_vars: # assume that if we have u we have v
    plot_uv = True
    if 'xi_u' in ds.u.dims:
        uv_grid = 'uv'
    elif 'xi_rho' in ds.u.dims:
        uv_grid = 'rho'

ndims = len(ds.salt.dims)
if ndims == 3:
    # like for a squeezed surface extraction
    s0 = ds.salt[nt,:,:].values
    if plot_uv:
        u0 = ds.u[nt,:,:].values
        v0 = ds.v[nt,:,:].values
elif ndims == 4:
    s0 = ds.salt[nt,-1,:,:].values
    if plot_uv:
        u0 = ds.u[nt,-1,:,:].values
        v0 = ds.v[nt,-1,:,:].values

# PLOTTING

plt.close('all')

pfun.start_plot(figsize=(10,10))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.pcolormesh(plon, plat, s0, cmap='Spectral_r')
pfun.add_coast(ax)
pfun.dar(ax)
pad = .02
ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])

if plot_uv:
    pfun.start_plot(figsize=(20,10))
    fig = plt.figure()

    ax = fig.add_subplot(121)
    if uv_grid == 'rho':
        ax.pcolormesh(plon, plat, u0, cmap='jet')
        ax.set_title('u rho-grid')
    elif uv_grid == 'uv':
        lon = ds.lon_u.values
        lat = ds.lat_u.values
        plon, plat = pfun.get_plon_plat(lon, lat)
        ax.pcolormesh(plon, plat, u0, cmap='jet')
        ax.set_title('u u-grid')
    pfun.add_coast(ax)
    pfun.dar(ax)
    pad = .02
    ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])

    ax = fig.add_subplot(122)
    if uv_grid == 'rho':
        ax.pcolormesh(plon, plat, v0, cmap='jet')
        ax.set_title('v rho-grid')
    elif uv_grid == 'uv':
        lon = ds.lon_v.values
        lat = ds.lat_v.values
        plon, plat = pfun.get_plon_plat(lon, lat)
        ax.pcolormesh(plon, plat, v0, cmap='jet')
        ax.set_title('v v-grid')
    pfun.add_coast(ax)
    pfun.dar(ax)
    pad = .02
    ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])

plt.show()
pfun.end_plot()

ds.close()