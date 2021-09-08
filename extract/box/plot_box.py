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
import matplotlib.pyplot as plt
import plotting_functions as pfun
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

plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)

s0 = ds.salt[nt,-1,:,:].values
rmask = ~np.isnan(s0) # True on water
plot_uv = False
if 'u' in ds.data_vars and 'v' in ds.data_vars:
    u0 = ds.u[nt,-1,:,:].values
    v0 = ds.v[nt,-1,:,:].values
    umask = ~np.isnan(u0)
    vmask = ~np.isnan(v0)
    plot_uv = True

# PLOTTING

plt.close('all')

if plot_uv:
    # show exact gridpoints and mask
    pfun.start_plot(figsize=(10,10))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    alpha=.2
    ms = 3
    ax.plot(ds.lon_rho.values, ds.lat_rho.values,'ok',alpha=alpha, ms=ms)
    ax.plot(ds.lon_u.values, ds.lat_u.values,'>g',alpha=alpha, ms=ms)
    ax.plot(ds.lon_v.values, ds.lat_v.values,'^r',alpha=alpha, ms=ms)
    ax.plot(ds.lon_rho.values[rmask], ds.lat_rho.values[rmask],'ok', ms=ms)
    ax.plot(ds.lon_u.values[umask], ds.lat_u.values[umask],'>g', ms=ms)
    ax.plot(ds.lon_v.values[vmask], ds.lat_v.values[vmask],'^r', ms=ms)
    pfun.add_coast(ax, color='b', linewidth=2)
    pfun.dar(ax)
    pad = .02
    ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])
    
    plt.show()
    pfun.end_plot()

if plot_uv:
    # quiver plot of velocity
    pfun.start_plot(figsize=(10,10))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # get velocity
    u0 = ds.u[nt,-1,:,:].values
    v0 = ds.v[nt,-1,:,:].values
    umask = ~np.isnan(u0)
    vmask = ~np.isnan(v0)
    # make interpolated u and v
    uu = u0.copy(); vv = v0.copy()
    uu[np.isnan(uu)] = 0
    vv[np.isnan(vv)] = 0
    UU = (uu[1:-1,1:]+uu[1:-1,:-1])/2
    VV = (vv[:1,1:-1] + vv[:-1,1:-1])/2
    UU[np.isnan(s0[1:-1,1:-1])] = np.nan
    VV[np.isnan(s0[1:-1,1:-1])] = np.nan

    alpha=.2
    ax.pcolormesh(plon, plat, s0, vmin=30, vmax=32, cmap='Spectral_r')
    ax.quiver(ds.lon_rho.values[1:-1,1:-1], ds.lat_rho.values[1:-1,1:-1], UU, VV)
    pfun.add_coast(ax, color='b', linewidth=2)
    pfun.dar(ax)
    pad = .02
    ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])
    
    plt.show()
    pfun.end_plot()

vn = 'oxygen'
if vn in ds.data_vars:
    # plot another variable
    vn = 'oxygen'
    if vn in ds.data_vars:
        v0 = ds[vn][nt,-1,:,:].values
    else:
        print('Variable %s not found' % (vn))
        sys.exit()
        
    if vn == 'oxygen':
        v0 = v0 * 32/1000
        vmin = 0
        vmax = 10
        cmap = cm.oxy
        tstr = 'Bottom DO (mg/L)'
    
    pfun.start_plot(figsize=(10,13))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    alpha=.2
    cs = ax.pcolormesh(plon, plat, v0, vmin=vmin, vmax=vmax, cmap=cmap)
    fig.colorbar(cs, ax=ax)
    pfun.add_coast(ax, color='k', linewidth=0.5)
    pfun.dar(ax)
    pfun.add_bathy_contours(ax, ds, depth_levs = [], txt=True)
    
    pad = .02
    ax.axis([plon[0,0]-pad, plon[-1,-1]+pad, plat[0,0]-pad, plat[-1,-1]+pad])
    
    ax.text(.05, .05, ot_dt[nt].strftime(tstr + '\n%Y.%m.%d\n%H:%M:%S' + ' UTC'), transform=ax.transAxes,
        bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        
    ax.set_title(box_name)
    
    plt.show()
    pfun.end_plot()

ds.close()