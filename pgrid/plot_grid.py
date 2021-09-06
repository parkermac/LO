# -*- coding: utf-8 -*-
"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.pu
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',
        type=str)
args = parser.parse_args()

from importlib import reload
import gfun
reload(gfun)
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

# **** USER CHOICES ********

# show plot of markers for rho, u, v grids with masking
flag_show_grids = False

# show resolution
flag_show_res = True

# show bathymetry on some sections
flag_show_sections = False

# **** END USER CHOICES ****

# select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] / fn
in_fn0 = Gr['gdir'] / 'grid_m00_r00_s00_x00.nc'

# load the data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
ds0 = xr.open_dataset(in_fn0)
z0 = -ds0.h.values
mask_rho = ds.mask_rho.values

plon, plat = pfun.get_plon(ds.lon_rho.values, ds.lat_rho.values)
buff = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-buff, plon[0,-1]+buff, plat[0,0]-buff, plat[-1,0]+buff)

zm = np.ma.masked_where(mask_rho == 0, z)
z0m = np.ma.masked_where(mask_rho == 0, z0)

# plotting
plt.close('all')

# set number of columns for plot 
NC = 1 # first guess
if flag_show_grids:
    NC += 1
    icg = NC       
if flag_show_res:
    NC += 1
    icr = NC
# this has to come last in the lineup
if flag_show_sections:
    NC += 1
    ics = NC
fig = plt.figure(figsize=(8*NC,8))

ax1 = fig.add_subplot(1,NC,1)
cmap1 = plt.get_cmap(name='gist_earth') # terrain, viridis
cs = ax1.pcolormesh(plon, plat, zm,
                   vmin=-300, vmax=10, cmap = cmap1)
fig.colorbar(cs, ax=ax1, extend='both')
pfun.add_coast(ax1)
pfun.dar(ax1)
ax1.axis(ax_lims)
ax1.set_title(Gr['gridname'] + '/' + fn)
ax1.text(.95, .05, str(mask_rho.shape), horizontalalignment='right',
         transform=ax1.transAxes)                   
#gfp.add_river_tracks(Gr, ds, ax1)
   
if flag_show_sections:
    color_list = ['orange', 'gold', 'greenyellow', 'lightgreen',
                    'aquamarine', 'cadetblue', 'royalblue', 'purple']
    lon_rho = ds['lon_rho'][:]
    lat_rho = ds['lat_rho'][:]
    NS = 8 # number of sections
    for ss in range(NS):
        ax = fig.add_subplot(NS,NC,ics*NS - ics*(ss+1) + ics)
        x = lon_rho[0, :]
        jj = int(lon_rho.shape[0]/(NS+1) * (ss+1))
        y = z[jj, :]/100
        y0 = z0[jj, :]/100
        ax.plot(x, y, '-r', linewidth=1)
        ax.plot(x, y0, '-c', linewidth=1)
        ax.plot(x, 0*x, '-', color=color_list[ss], linewidth=1)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(-5, 1)
        if ss == NS-1:
            ax.set_title('Z/(100 m)')
        
        ax1.plot([x[0], x[-1]], [lat_rho[jj,0], lat_rho[jj, -1]], '-',
            color=color_list[ss], linewidth=2)
        
if flag_show_grids:
    # NOTE: you need to have run make_extras.py for this to work
    lon_dict = dict()
    lat_dict = dict()
    mask_dict = dict()
    tag_list = ['rho', 'u', 'v', 'psi']
    for tag in tag_list:
        lon_dict[tag] = ds.variables['lon_'+tag][:]
        lat_dict[tag] = ds.variables['lat_'+tag][:]
        mask_dict[tag] = ds.variables['mask_'+tag][:]
    marker_dict = {'rho': 'ok',
                 'u': '>r',
                 'v': '^b',
                 'psi': 'xg'}
    
    ax = fig.add_subplot(1,NC,icg)
    for tag in tag_list:
        ax.plot(lon_dict[tag][mask_dict[tag]==1], lat_dict[tag][mask_dict[tag]==1],
                marker_dict[tag])
        ax.plot(lon_dict[tag][mask_dict[tag]==0], lat_dict[tag][mask_dict[tag]==0],
                marker_dict[tag], alpha = .2)
    pfun.dar(ax)
    ax.set_xlim(ax_lims[:2])
    ax.set_ylim(ax_lims[-2:])
    
    pfun.add_coast(ax)
    ax.axis(ax_lims)
    gfp.add_river_tracks(Gr, ds, ax)
    
if flag_show_res:
    DX = 1/ds['pm'][:]
    DY = 1/ds['pn'][:]
    res = np.maximum(DX, DY)
    ax = fig.add_subplot(1,NC,icr)
    cmap = plt.get_cmap(name='rainbow_r') # terrain, viridis
    cs = ax.pcolormesh(plon, plat, res,
                       vmin=500, vmax=1500, cmap = cmap)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(ax_lims)
    ax.set_title('Grid Resolution - maximum (m)')

ds.close()
ds0.close()

plt.show()
