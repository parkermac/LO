# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:45:54 2016

@author: PM5

Smooth the grid.

"""

from importlib import reload
import gfun
reload(gfun)
Gr =gfun.gstart()
import gfun_utility as gfu
reload(gfu)
import pfun

import netCDF4 as nc
import os
import shutil
import numpy as np
import pickle

# load the default choices
dch = pickle.load(open(Gr['gdir'] + 'choices.p', 'rb'))

#%% select grid file
fn = gfun.select_file(Gr)
in_fn = Gr['gdir'] + fn
# create new file name
fn_new = gfun.increment_filename(fn, tag='_s')
out_fn = Gr['gdir'] + fn_new

# retrieve the output

ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]
plon = ds.variables['lon_psi_ex'][:]
plat = ds.variables['lat_psi_ex'][:]
lon = ds.variables['lon_rho'][:]
lat = ds.variables['lat_rho'][:]
dx = 1/ds.variables['pm'][:]
dy = 1/ds.variables['pn'][:]
ds.close()

#%% Smoothing set up

MSK = mask_rho
Hobs = -z.copy()

rx0max = 0.15

# Make sure that anything NOT MASKED is not shallower than min_depth.
if dch['use_min_depth']:
    Hobs[(MSK==1) & (Hobs < dch['min_depth'])] = dch['min_depth']

# create the area matrix
AreaMatrix = dx * dy

# create smoothed bathymetry

import time
tt0 = time.time()

rx0max = 0.15

# Shift the whole depth grid for the case where we have active bathymetry
# on land (e.g. when using wet_dry), becasue otherwise the smoothing
# code will fail.  We shift it back at the end.  All the shifting is done in
# gfu.GRID_PlusMinusScheme_rx0().
if dch['min_depth'] > 0:
    shift = 0
elif dch['min_depth'] <= 0:
    shift = dch['min_depth'] - 1
    
# Do the smoothing.
Hnew = gfu.GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix,
            fjord_cliff_edges=True, shift=shift)

print('Smoothing took %0.1f seconds' % (time.time() - tt0))

# Again, make sure that anything NOT MASKED is not shallower than min_depth.
if dch['use_min_depth']:
    Hnew[(MSK==1) & (Hnew < dch['min_depth'])] = dch['min_depth']

zn = -Hnew

# Save the output file

print('Creating ' + out_fn)
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist
shutil.copyfile(in_fn, out_fn)
ds = nc.Dataset(out_fn, 'a')
    
ds['h'][:] = -zn
ds['mask_rho'][:] = MSK

#%% plotting

if True:

    import matplotlib.pyplot as plt

    ax_lims = (plon[0,0], plon[0,-1], plat[0,0], plat[-1,0])

    z_m = np.ma.masked_where(mask_rho==0, z)
    zn_m = np.ma.masked_where(mask_rho==0, zn)
    dz = zn - z

    plt.close()

    fig = plt.figure(figsize=(14,10))

    ax = fig.add_subplot(121)
    cmap1 = plt.get_cmap(name='terrain')
    cs = ax.pcolormesh(plon, plat, zn_m,
                       vmin=-200, vmax=200, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title(Gr['gridname'] + '/' + fn_new)

    ax = fig.add_subplot(122)
    cmap1 = plt.get_cmap(name='bwr')
    cs = ax.pcolormesh(plon, plat, dz,
                       vmin=-10, vmax=10, cmap = cmap1)
    fig.colorbar(cs, ax=ax, extend='both')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(pfun.get_aa(ds))
    ax.set_title('New Z - Old Z (m)')

    plt.show()

ds.close()
